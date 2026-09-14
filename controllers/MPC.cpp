class MPC {
public:
    struct State {
        double y;    // lateral error (e_y) [m]
        double v_y;  // lateral velocity [m/s]
        double psi;  // heading error (e_psi) [rad]
        double r;    // yaw rate [rad/s]
    };

    struct Control {
        double delta; // steering angle [rad]
    };

    struct Params {
        double m;    // mass [kg]
        double I_z;  // yaw inertia [kg·m^2]
        double l_f;  // distance CoG -> front axle [m]
        double l_r;  // distance CoG -> rear axle [m]
        double C_f;  // front cornering stiffness in negative [N/rad]
        double C_r;  // rear cornering stiffness in negative [N/rad]
    };

     static void buildLinearModel(const State& x, const Control& u, 
                                 Eigen::Matrix<double, 4, 4>& A,
                                 Eigen::Matrix<double, 4, 1>& B,
                                 Eigen::Matrix<double, 4, 1>& c,
                                 const Params& p, double v_x, double kappa_ref)
    {
        // --- Parse state and input ---
        double y     = x.y;
        double v_y   = x.v_y;
        double psi   = x.psi;
        double r     = x.r;
        double delta = u.delta;

        // --- Initialize matrices ---
        A.setZero();
        B.setZero();
        c.setZero();

        // --- A matrix ---
        A(0,0) = 0.0;
        A(0,1) = 1.0;
        A(0,2) = v_x;
        A(0,3) = 0.0;

        A(1,0) = 0.0;
        A(1,1) = (p.C_f + p.C_r) / (p.m * v_x);
        A(1,2) = 0.0;
        A(1,3) = (p.l_f * p.C_f - p.l_r * p.C_r) / (p.m * v_x) - v_x;

        A(2,0) = 0.0;
        A(2,1) = 0.0;
        A(2,2) = 0.0;
        A(2,3) = 1.0;

        A(3,0) = 0.0;
        A(3,1) = (p.l_f * p.C_f - p.l_r * p.C_r) / (p.I_z * v_x);
        A(3,2) = 0.0;
        A(3,3) = (p.l_f * p.l_f * p.C_f + p.l_r * p.l_r * p.C_r) / (p.I_z * v_x);


        // --- B matrix ---
        B(0,0) = 0.0;
        B(1,0) = -p.C_f / p.m;
        B(2,0) = 0.0;
        B(3,0) = -p.l_f * p.C_f / p.I_z;

        //xdot (DYNAMIC BICCYLE MODEL EQUATIONS BELOW)
        State f;

        // Slip angles
        double alpha_f = -u.delta + std::atan2((v_y + p.l_f * r), v_x);
        double alpha_r = std::atan2((v_y - p.l_r * r), v_x);

        // Tire forces
        double Fy_f = p.C_f * alpha_f;
        double Fy_r = p.C_r * alpha_r;

        // Yaw moment
        double Mz = p.l_f * Fy_f * std::cos(u.delta) - p.l_r * Fy_r;

        // Continuous dynamics (xdot)
        f.y   = v_y;
        f.v_y = (Fy_f * std::cos(u.delta) + Fy_r) / p.m - v_x * r;
        // Heading-error dynamics in curvilinear coordinates.
        f.psi = r - v_x * kappa_ref;
        f.r   = Mz / p.I_z;

        Eigen::Matrix<double, 4, 1> x_vec;
        x_vec << y, v_y, psi, r;

        Eigen::Matrix<double, 1, 1> u_vec;
        u_vec << delta;

        Eigen::Matrix<double, 4, 1> f_vec;
        f_vec << f.y, f.v_y, f.psi, f.r;

        c = f_vec - A * x_vec - B * u_vec;
    }

    static std::tuple<Eigen::VectorXd, std::vector<double>, std::vector<std::pair<double,double>>>
    buildReference(const utfr_msgs::msg::ParametricSpline& spline, double dt, int N, int NX, double v_x, const utfr_msgs::msg::EgoState &ego_state_new, const utfr_msgs::msg::EgoState &ego_at_path_time, double &s0){

        // Ego offset in GLOBAL frame
        double dx_global = ego_state_new.pose.pose.position.x - ego_at_path_time.pose.pose.position.x;
        double dy_global = ego_state_new.pose.pose.position.y - ego_at_path_time.pose.pose.position.y;

        // Use existing util helper to get yaw from quaternion
        double yaw_new = util::quaternionToYaw(ego_state_new.pose.pose.orientation);
        double yaw_ref = util::quaternionToYaw(ego_at_path_time.pose.pose.orientation);
        double dpsi = std::atan2(std::sin(yaw_new - yaw_ref), std::cos(yaw_new - yaw_ref)); // Angle Wrapping

        // Convert offset from GLOBAL frame to ego_at_path_time LOCAL frame
        double cos_ref = std::cos(yaw_ref);
        double sin_ref = std::sin(yaw_ref);
        double dx_local = cos_ref * dx_global + sin_ref * dy_global;
        double dy_local = -sin_ref * dx_global + cos_ref * dy_global;

        const std::vector<double>& x_params = spline.x_params;
        const std::vector<double>& y_params = spline.y_params;

        const double dt_pub = 0.01; // time step for updating s0 based on current velocity and MPC publish rate

        // === Reference containers ===
        Eigen::VectorXd x_ref_vec(NX * N);
        std::vector<double> delta_ref(N, 0.0);
        std::vector<std::pair<double, double>> path_points; //for viz
        path_points.reserve(N);

        // === Vehicle geometry ===
        const double Lf = 0.811;  // front axle to CoG [m]
        const double Lr = 0.719;  // rear axle to CoG [m]
        const double L  = Lf + Lr; // wheelbase

        // === computes the new starting point for the available spline ===
        {
          // Horner's Method for 1st derivative
          double dx_ds0 = (((5*x_params[0]*s0 + 4*x_params[1])*s0 + 3*x_params[2])*s0 + 2*x_params[3])*s0 + x_params[4];
          double dy_ds0 = (((5*y_params[0]*s0 + 4*y_params[1])*s0 + 3*y_params[2])*s0 + 2*y_params[3])*s0 + y_params[4];

          double deriv_mag0 = std::sqrt(dx_ds0*dx_ds0 + dy_ds0*dy_ds0);
          
          double ds_tick = 0.0;

          if (deriv_mag0 > 1e-6) {
              ds_tick = (std::max(0.0, v_x) * dt_pub) / deriv_mag0; // normalize
          }

          double max_ds_tick = 1/N;
          ds_tick = std::clamp(ds_tick, 0.0, max_ds_tick); // safety
          s0 = std::clamp(s0 + ds_tick, 0.0, 1.0);
        }

        double s = s0;
        // === from s0, discretize along the spline and compute reference states and inputs ===
        for (int i = 0; i < N; ++i){
              // --- Compute local derivatives wrt s ---
              // Horner's Method
              double dx_ds = (((5*x_params[0]*s + 4*x_params[1])*s + 3*x_params[2])*s + 2*x_params[3])*s + x_params[4];
              double dy_ds = (((5*y_params[0]*s + 4*y_params[1])*s + 3*y_params[2])*s + 2*y_params[3])*s + y_params[4];

              // --- Compute local arc-length scaling ---
              double d2 = dx_ds*dx_ds + dy_ds*dy_ds;
              double deriv_mag = std::sqrt(d2);

              double max_ds = 0.04; //max step along spline to avoid overshoot
              double ds_arc = (v_x * dt) / deriv_mag;
              ds_arc = std::clamp(ds_arc, 0.0, max_ds);

              // Distance of x and y params from car
              // Horner's Method
              double x_path = (((((x_params[0]*s + x_params[1])*s + x_params[2])*s + x_params[3])*s + x_params[4])*s + x_params[5]);
              double y_path = (((((y_params[0]*s + y_params[1])*s + y_params[2])*s + y_params[3])*s + y_params[4])*s + y_params[5]);

              // Spline point (x_path, y_path) is in ego_at_path_time local frame
              // Subtract offset: car moved forward, so path appears to move backward relative to car
              double x_shifted_local = x_path - dx_local;
              double y_shifted_local = y_path - dy_local;

              // Rotate from ego_at_path_time local to ego_state_new local by dpsi
              double cos_dpsi = std::cos(dpsi);
              double sin_dpsi = std::sin(dpsi);
              double x_shifted = cos_dpsi * x_shifted_local - sin_dpsi * y_shifted_local;
              double y_shifted = sin_dpsi * x_shifted_local + cos_dpsi * y_shifted_local;

              path_points.emplace_back(x_shifted, y_shifted);

              // Horner's Method
              double ddx_ds = ((20*x_params[0]*s + 12*x_params[1])*s + 6*x_params[2])*s + 2*x_params[3];
              double ddy_ds = ((20*y_params[0]*s + 12*y_params[1])*s + 6*y_params[2])*s + 2*y_params[3];

              // Curvature: Optimized to use existing d2 and deriv_mag
              double denom = d2 * deriv_mag; 
              double kappa = (denom > 1e-10) ? (dx_ds * ddy_ds - dy_ds * ddx_ds) / denom : 0;

              double max_kappa = std::tan(1.5) / L;
              kappa = std::clamp(kappa, -max_kappa, max_kappa);

              // === Reference steering from curvature ===
              delta_ref[i] = std::atan(L * kappa);

              // === Reference yaw rate ===
              double r_ref = v_x * kappa;
              double max_r = 5.0;  // rad/s absolute safety limit
              r_ref = std::clamp(r_ref, -max_r, max_r);

              // === Pack state reference vector [e_y, v_y, e_psi, r] ===
              // UPDATED: Tracking zero error
              x_ref_vec.segment(i * NX, NX) << 0.0, 0.0, 0.0, r_ref;

              // --- Advance along spline ---
              s = std::clamp(s + ds_arc, 0.0, 1.0);
          }

          return std::make_tuple(x_ref_vec, delta_ref, path_points);
      }
};

void ControllerNode::setPrevD(double d){
    delta_prev_ = d;
}

double ControllerNode::getPrevD(){
    return delta_prev_;
}

utfr_msgs::msg::TargetState
ControllerNode::LTVMPC(const utfr_msgs::msg::ParametricSpline& spline, utfr_msgs::msg::VelocityProfile &velocity_profile, utfr_msgs::msg::EgoState &ego_state_new){

    // === Extract vehicle velocity from ego_state_new === (safe guard)
    double v_x = ego_state_new.vel.twist.linear.x;
    // minimum velocity to avoid numerical issues in linearization and QP
    if (v_x < 3.5){
      v_x = 3.5;
    }

    // === Calculate Initial Errors at s0 ===
    double dx_ds0 = (((5*spline.x_params[0]*s0 + 4*spline.x_params[1])*s0 + 3*spline.x_params[2])*s0 + 2*spline.x_params[3])*s0 + spline.x_params[4];
    double dy_ds0 = (((5*spline.y_params[0]*s0 + 4*spline.y_params[1])*s0 + 3*spline.y_params[2])*s0 + 2*spline.y_params[3])*s0 + spline.y_params[4];

    double x_path_s0 = (((((spline.x_params[0]*s0 + spline.x_params[1])*s0 + spline.x_params[2])*s0 + spline.x_params[3])*s0 + spline.x_params[4])*s0 + spline.x_params[5]);
    double y_path_s0 = (((((spline.y_params[0]*s0 + spline.y_params[1])*s0 + spline.y_params[2])*s0 + spline.y_params[3])*s0 + spline.y_params[4])*s0 + spline.y_params[5]);

    // Convert spline starting point to car's CURRENT local frame
    double yaw_new = util::quaternionToYaw(ego_state_new.pose.pose.orientation);
    double yaw_ref = util::quaternionToYaw(ego_at_path_time_.pose.pose.orientation);
    double dpsi = std::atan2(std::sin(yaw_new - yaw_ref), std::cos(yaw_new - yaw_ref));

    double dx_global = ego_state_new.pose.pose.position.x - ego_at_path_time_.pose.pose.position.x;
    double dy_global = ego_state_new.pose.pose.position.y - ego_at_path_time_.pose.pose.position.y;
    double cos_ref = std::cos(yaw_ref); double sin_ref = std::sin(yaw_ref);
    double dx_local = cos_ref * dx_global + sin_ref * dy_global;
    double dy_local = -sin_ref * dx_global + cos_ref * dy_global;

    double x_shifted_local = x_path_s0 - dx_local;
    double y_shifted_local = y_path_s0 - dy_local;

    double cos_dpsi = std::cos(dpsi); double sin_dpsi = std::sin(dpsi);
    double y_err_initial = sin_dpsi * x_shifted_local + cos_dpsi * y_shifted_local; 

    double psi_path_s0 = std::atan2(dy_ds0, dx_ds0) + dpsi;
    double psi_err_initial = -psi_path_s0; // Car is at yaw 0 in its own frame

    // === Vehicle parameters ===
    MPC::Params params;
    params.m   = 200.0;                      // sprung_mass [kg]
    params.I_z = 110.0;                      // Izz [kg*m^2]
    params.l_f = 0.811;                      // a_cg - distance from CG to front axle [m]
    params.l_r = 0.719;                      // b_cg - distance from CG to rear axle [m]
    params.C_f = -560.0 * (180.0 / M_PI);    // C_f tire cornering coefficient [N/rad]
    params.C_r = -560.0 * (180.0 / M_PI);    // C_r tire cornering coefficient [N/rad]

    // === Initial state (Error formulation) ===
    MPC::State x0;
    x0.y   = -y_err_initial; // e_y
    x0.v_y = ego_state_new.vel.twist.linear.y;
    x0.psi = psi_err_initial; // e_psi
    // initial yaw rate from prev. steering.
    double len = params.l_f + params.l_r;
    x0.r   = (v_x * getPrevD()) / ((len) + (((params.m * v_x * v_x) / len)) * ((params.l_f/params.C_f)-(params.l_r/params.C_r))); 

    // === Allocate QP matrices ===
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nVar);
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(nConstraintsX, nVar);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(nConstraintsX);
    Eigen::MatrixXd Aineq = Eigen::MatrixXd::Zero(nConstraintsU, nVar);
    Eigen::VectorXd lbAineq = Eigen::VectorXd::Zero(nConstraintsU);
    Eigen::VectorXd ubAineq = Eigen::VectorXd::Zero(nConstraintsU);

    // Anchor the first predicted state to the measured initial error state.
    Aeq.block(0, 0, NX, NX) = Eigen::MatrixXd::Identity(NX, NX);
    beq.segment(0, NX) << x0.y, x0.v_y, x0.psi, x0.r;

    // === Build reference trajectory ===
    Eigen::VectorXd x_ref_vec;
    std::vector<double> delta_ref;
    std::vector<std::pair<double, double>> path_points;
    std::tie(x_ref_vec, delta_ref, path_points) = MPC::buildReference(spline, dt_, N, NX, v_x, ego_state_new, ego_at_path_time_, s0);

    // Build gradient vector g
    g.setZero(nVar);
    g.head(N * NX) = -2.0 * Qblk * x_ref_vec;

    //Previous control input
    double delta_prev = getPrevD();
    Eigen::VectorXd d_prev = Eigen::VectorXd::Zero(N);
    d_prev(0) = delta_prev;
    g.segment(N * NX, N * NU) = -2.0 * D.transpose() * Rblk * d_prev;

    // === Visualize reference and actual points (from prev. iteration) === 
    visualization_msgs::msg::Marker ref_points_marker;
    visualization_msgs::msg::Marker prev_mpc_marker;
    // convert frenet frame into baselink
    std::vector<std::pair<double, double>> prev_mpc_points;
    if (has_prev_mpc_solution_) {
      prev_mpc_points.reserve(N);
      auto pathHeading = [&](int idx) {
        const int prev_idx = std::max(0, idx - 1);
        const int next_idx = std::min(N - 1, idx + 1);
        const double dx = path_points[next_idx].first - path_points[prev_idx].first;
        const double dy = path_points[next_idx].second - path_points[prev_idx].second;
        return std::atan2(dy, dx);
      };
      for (int k = 0; k < N; ++k) {
        const double e_y = x_opt_[k * NX + 0];
        const double yaw_path = pathHeading(k);
        const double normal_x = -std::sin(yaw_path);
        const double normal_y = std::cos(yaw_path);

        const double x = path_points[k].first + e_y * normal_x;
        const double y = path_points[k].second + e_y * normal_y;
        prev_mpc_points.emplace_back(x, y);
      }
    }
    vizPoints(ref_points_marker, prev_mpc_marker, N, path_points, prev_mpc_points);

    // === Linearization around reference ===
    for (int i = 0; i < N-1; ++i) {
        Eigen::Matrix<double, NX, NX> A;
        Eigen::Matrix<double, NX, NU> B;
        Eigen::Matrix<double, NX, 1> c;

        // Linearize at reference state and input
        Eigen::VectorXd x_ref_i = x_ref_vec.segment(i*NX, NX);
        MPC::State xbar_ref{ x_ref_i(0), x_ref_i(1), x_ref_i(2), x_ref_i(3) };
        MPC::Control ubar_ref{ delta_ref[i] };

        // Extract kappa (since r_ref = v_x * kappa)
        double kappa_ref = x_ref_i(3) / v_x;

        MPC::buildLinearModel(xbar_ref, ubar_ref, A, B, c, params, v_x, kappa_ref);

        // Discretize using dt
        Eigen::Matrix<double, NX, NX> A_d = Eigen::Matrix<double, NX, NX>::Identity() + A*dt_;
        Eigen::Matrix<double, NX, NU> B_d = B*dt_;
        Eigen::Matrix<double, NX, 1> c_d = c*dt_;

        // Fill equality constraints: x_{k+1} = A x_k + B u_k + c
        const int row = (i + 1) * NX;
        Aeq.block(row, i*NX, NX, NX)       = -A_d;
        Aeq.block(row, (i+1)*NX, NX, NX)   = Eigen::MatrixXd::Identity(NX, NX);
        Aeq.block(row, N*NX + i*NU, NX, NU) = -B_d;
        beq.segment(row, NX) = c_d;
    }
    // Fill inequality constraints for steering rate limits: |delta_k - delta_{k-1}| <= max_delta_rate
    //      - assumes we can go from left lock to right lock in 0.5s
    Aineq.block(0, N*NX, nConstraintsU, N*NU) = D;
    lbAineq = Eigen::VectorXd::Constant(nConstraintsU, -dt_*1.4); // left -> right takes 1s
    lbAineq(0) = -max_steering_angle_; // rely on PI controllers in sim instead of bounding initial steering input to a range.
    ubAineq = Eigen::VectorXd::Constant(nConstraintsU, dt_*1.4);
    ubAineq(0) =  max_steering_angle_; 

    // === Variable bounds (steering) ===
    Eigen::VectorXd lb = -1e5 * Eigen::VectorXd::Ones(nVar);
    Eigen::VectorXd ub =  1e5 * Eigen::VectorXd::Ones(nVar);
    for (int i = 0; i < N; ++i) {
        int idu = N * NX + i * NU;
        lb(idu) = -max_steering_angle_; // min steering [rad]
        ub(idu) =  max_steering_angle_; // max steering [rad]
        int idx = i * NX;
        double len = params.l_f + params.l_r;
        double max_yaw_rate = (v_x * max_steering_angle_) / ((len) + (((params.m * v_x * v_x) / len)) * ((params.l_f/params.C_f)-(params.l_r/params.C_r))); 
        lb(idx + 3) = -max_yaw_rate; // max yaw rate  [rad/s] (velocity dependent max yaw rate)
        ub(idx + 3) =  max_yaw_rate; // min yaw rate  [rad/s]
    }

    // === Solve QP ===
    if(!has_prev_mpc_solution_) {
      proxqp_->init(H, g, Aeq, beq, Aineq, lbAineq, ubAineq); // start new solution
      proxqp_->solve();
    } else {
      proxqp_->update(std::nullopt, g, Aeq, beq, Aineq, lbAineq, ubAineq); // warm start with previous solution
      proxqp_->solve();
    }
    auto result = proxqp_->results.info.status;
    //Initialize TargetState
    utfr_msgs::msg::TargetState target;
    if (result == proxsuite::proxqp::QPSolverOutput::PROXQP_SOLVED) {
        Eigen::Map<Eigen::VectorXd>(x_opt_, nVar) = proxqp_->results.x;        
        has_prev_mpc_solution_ = true;
        double steering = x_opt_[N * NX]; // first control input (steering)
        steering = std::clamp(steering, -max_steering_angle_, max_steering_angle_); //safety clamp max steering
        target.steering_angle = steering;
        setPrevD(steering);
        RCLCPP_INFO(this->get_logger(), "Optimal steering = %f. solve time: %f", 
                    steering, proxqp_->results.info.run_time / 1000.0); // solve time in ns
    } else {
        RCLCPP_WARN(this->get_logger(), "QP failed: %d", result);
        has_prev_mpc_solution_ = false;
        // Get primal solution even on failure
        Eigen::VectorXd z = proxqp_->results.x;
        if(solver_debug_) {
            RCLCPP_ERROR(this->get_logger(), "QP solver failed with code %d", result);

            Eigen::VectorXd Ax_eq = Aeq * z;
            Eigen::VectorXd Ax_ineq = Aineq * z;
            RCLCPP_ERROR(this->get_logger(), "=== CONSTRAINT VIOLATION REPORT ===");
            for (int i = 0; i < nConstraintsX; ++i) {
                double viol = std::abs(Ax_eq(i) - beq(i));
                if (viol > 1e-4) {
                    int step = i / NX;
                    int state = i % NX;
                    const char* snames[] = {"e_y","v_y","e_psi","r"};
                    RCLCPP_ERROR(this->get_logger(),
                        "  EQ[step=%d, %s]: Ax=%.4f, beq=%.4f, viol=%.4f",
                        step, snames[state], Ax_eq(i), beq(i), viol);
                }
            }
            for (int i = 0; i < nConstraintsU; ++i) {
                double lb_viol = lbAineq(i) - Ax_ineq(i);
                double ub_viol = Ax_ineq(i) - ubAineq(i);
                if (lb_viol > 1e-4 || ub_viol > 1e-4) {
                    RCLCPP_ERROR(this->get_logger(),
                        "  INEQ[u_row=%d]: Ax=%.4f, lb=%.4f, ub=%.4f, lb_viol=%.4f, ub_viol=%.4f",
                        i, Ax_ineq(i), lbAineq(i), ubAineq(i), lb_viol, ub_viol);
                }
            }
          }
    }

    splitVelocity(target, velocity_profile, ego_state_new);

    return target;

}
