class MPC {
public:
    struct State {
        double y;    // lateral displacement [m]
        double v_y;  // lateral velocity [m/s]
        double psi;  // yaw angle [rad]
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
                                 const Params& p, double v_x)
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
        f.psi = r;
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
    buildReference(const spline_msgs::msg::ParametricSpline& spline, int N, int NX, double v_x){

        std::vector<double> x_params = spline.x_params;
        std::vector<double> y_params = spline.y_params;

        // === Fixed timestep for MPC ===
        double dt = 0.04;

        // === Reference containers ===
        Eigen::VectorXd x_ref_vec(NX * N);
        std::vector<double> delta_ref(N, 0.0);
        std::vector<std::pair<double, double>> path_points; //for viz
        path_points.reserve(N);

        // === Vehicle geometry ===
        const double Lf = 0.811;  // front axle to CoG [m]
        const double Lr = 0.719;  // rear axle to CoG [m]
        const double L  = Lf + Lr; // wheelbase

        double s = 0.0;

        for (int i = 0; i < N; ++i){
              // --- Compute local derivatives wrt s ---
              double dx_ds = 5*x_params[0]*pow(s,4) + 4*x_params[1]*pow(s,3) +
                            3*x_params[2]*pow(s,2) + 2*x_params[3]*s + x_params[4];
              double dy_ds = 5*y_params[0]*pow(s,4) + 4*y_params[1]*pow(s,3) +
                            3*y_params[2]*pow(s,2) + 2*y_params[3]*s + y_params[4];

              // --- Compute local arc-length scaling ---
              double deriv_mag = std::sqrt(dx_ds*dx_ds + dy_ds*dy_ds);

              double max_ds = 0.04; //max step along spline to avoid overshoot
              double ds_arc = (v_x * dt) / deriv_mag;
              ds_arc = std::clamp(ds_arc, 0.0, max_ds);

              //Distance of x and y params from car
              double x_path = x_params[0] * pow(s, 5) + x_params[1] * pow(s, 4) +
                        x_params[2] * pow(s, 3) + x_params[3] * pow(s, 2) +
                        x_params[4] * s + x_params[5];

              double y_path = y_params[0] * pow(s, 5) + y_params[1] * pow(s, 4) +
                        y_params[2] * pow(s, 3) + y_params[3] * pow(s, 2) +
                        y_params[4] * s + y_params[5];

              path_points.emplace_back(x_path, y_path);

              double ddx_ds = 20*x_params[0]*pow(s,3) + 12*x_params[1]*pow(s,2) + 6*x_params[2]*s + 2*x_params[3];
              double ddy_ds = 20*y_params[0]*pow(s,3) + 12*y_params[1]*pow(s,2) + 6*y_params[2]*s + 2*y_params[3];


              // === Compute path heading and curvature ===
              double psi_desired = atan2(dy_ds, dx_ds);
              
              double denom = std::pow(dx_ds * dx_ds + dy_ds * dy_ds, 1.5);
              double kappa = (denom > 1e-10) ? (dx_ds * ddy_ds - dy_ds * ddx_ds) / denom : 0;
              //kappa = std::clamp(kappa, -0.2, 0.2);

              double max_kappa = tan(0.7) / L;  //max angle of 0.34rads
              kappa = std::clamp(kappa, -max_kappa, max_kappa);


              // === Reference steering from curvature ===
              delta_ref[i] = std::atan(L * kappa);

              // === Reference yaw rate ===
              double r_ref = v_x * kappa;
              double max_r = 5.0;  // rad/s absolute safety limit
              r_ref = std::clamp(r_ref, -max_r, max_r);

              // === Reference lateral velocity ===
              double v_y_ref = 0;

              // === Pack state reference vector [y, v_y, psi, r] ===
              x_ref_vec.segment(i * NX, NX) << y_path, v_y_ref, psi_desired, r_ref;

              // --- Advance along spline ---
              s += ds_arc;
              s = std::clamp(s, 0.0, 1.0);
          }

          return std::make_tuple(x_ref_vec, delta_ref, path_points);
      }
};

utfr_msgs::msg::TargetState
ControllerNode::LTVMPC(const utfr_msgs::msg::ParametricSpline& spline, utfr_msgs::msg::VelocityProfile &velocity_profile) {

  // === Extract vehicle velocity === (safe guard)
    double v_x = ego_state_->vel.twist.linear.x;
    if (v_x < 1.0){
      v_x = 1.0;
    }

    // === MPC setup ===
    constexpr int N = 25;           // Horizon length
    constexpr int NX = 4;           // y, v_y, psi, r
    constexpr int NU = 1;           // Steering
    constexpr int nVar = N*NX + N*NU;        // Total variables (states + controls)
    constexpr int nConstraints = (N-1)*NX;   // Dynamics constraints

    // === Initial state (car's local frame) ===
    MPC::State x0;
    x0.y   = 0.0;  // Car at origin
    x0.v_y = ego_state_->vel.twist.linear.y;
    x0.psi = 0.0;
    x0.r   = ego_state_->vel.twist.angular.z;

    // === Vehicle parameters ===
    MPC::Params params;
    params.m   = 200.0;                      // sprung_mass [kg]
    params.I_z = 110.0;                      // Izz [kg*m^2]
    params.l_f = 0.811;                      // a_cg - distance from CG to front axle [m]
    params.l_r = 0.719;                      // b_cg - distance from CG to rear axle [m]
    params.C_f = -560.0 * (180.0 / M_PI);    // C_f tire cornering coefficient [N/rad]
    params.C_r = -560.0 * (180.0 / M_PI);    // C_r tire cornering coefficient [N/rad]

    // === Cost matrices ===
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(NX, NX);
    Q(0,0) = q_y_;  // lateral position
    Q(1,1) = q_vy_;   // lateral velocity
    Q(2,2) = q_psi_;  // heading
    Q(3,3) = q_r_;   // yaw rate

    // === Allocate QP matrices ===
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nVar, nVar);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nVar);
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(nConstraints, nVar);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(nConstraints);

    //Matrix Qblk contains matrix Q in diagonal
    Eigen::MatrixXd Qblk = Eigen::MatrixXd::Zero(N * NX, N * NX);
    for (int i = 0; i < N; ++i)
        Qblk.block(i * NX, i * NX, NX, NX) = Q;

    //D matrix (change in control input penalty) [u0, u1 - u0, u2 - u1]
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        D(i, i) = 1.0;
        if (i > 0)
            D(i, i - 1) = -1.0;
    }

    // --- Input weighting ---
    double R_val = r_delta_;
    Eigen::MatrixXd Rblk = R_val * Eigen::MatrixXd::Identity(N, N);

    //Build Hessian H (apply 2x scaling since qpOASES uses 0.5)
    H.setZero(nVar, nVar);
    H.block(0, 0, N * NX, N * NX) = 2.0 * Qblk;
    Eigen::MatrixXd H_u = 2.0 * (D.transpose() * Rblk * D);
    H.block(N * NX, N * NX, N * NU, N * NU) = H_u;

    // regularize and symmetrize (for numerical stability)
    H += 1e-8 * Eigen::MatrixXd::Identity(nVar, nVar);
    H = 0.5 * (H + H.transpose());

    // === Build reference trajectory ===
    Eigen::VectorXd x_ref_vec;
    std::vector<double> delta_ref;
    std::vector<std::pair<double, double>> path_points;
    std::tie(x_ref_vec, delta_ref, path_points) = MPC::buildReference(*path_, N, NX, v_x);

    // Log first few reference points
    RCLCPP_INFO(this->get_logger(),
        "Ref[0]: y=%.3f, vy=%.3f, psi=%.3f, r=%.3f, delta=%.3f",
        x_ref_vec(0), x_ref_vec(1), x_ref_vec(2), x_ref_vec(3), delta_ref[0]);

    RCLCPP_INFO(this->get_logger(),
        "Ref[1]: y=%.3f, vy=%.3f, psi=%.3f, r=%.3f, delta=%.3f",
        x_ref_vec(4), x_ref_vec(5), x_ref_vec(6), x_ref_vec(7), delta_ref[1]);

    // Build gradient vector g
    g.setZero(nVar);
    g.head(N * NX) = -2.0 * Qblk * x_ref_vec;

    //Prevoius control input
    double delta_prev = getPrevD();
    Eigen::VectorXd d_prev = Eigen::VectorXd::Zero(N);
    d_prev(0) = delta_prev;
    g.segment(N * NX, N * NU) = -2.0 * D.transpose() * Rblk * d_prev;

    double dt = 0.04;

    // Visualize all linearization reference points
    visualization_msgs::msg::Marker ref_points_marker;
    vizPoints(ref_points_marker, N, path_points);

    // === Linearization around reference ===
    for (int i = 0; i < N-1; ++i) {
        Eigen::Matrix<double, NX, NX> A;
        Eigen::Matrix<double, NX, NU> B;
        Eigen::Matrix<double, NX, 1> c;

        // Linearize at reference state and input
        Eigen::VectorXd x_ref_i = x_ref_vec.segment(i*NX, NX);
        MPC::State xbar_ref{ x_ref_i(0), x_ref_i(1), x_ref_i(2), x_ref_i(3) };
        MPC::Control ubar_ref{ delta_ref[i] };

        MPC::buildLinearModel(xbar_ref, ubar_ref, A, B, c, params, v_x);

        // Discretize using dt
        Eigen::Matrix<double, NX, NX> A_d = Eigen::Matrix<double, NX, NX>::Identity() + A*dt;
        Eigen::Matrix<double, NX, NU> B_d = B*dt;
        Eigen::Matrix<double, NX, 1> c_d = c*dt;

        // Fill equality constraints: x_{k+1} = A x_k + B u_k + c
        Aeq.block(i*NX, i*NX, NX, NX)       = -A_d;
        Aeq.block(i*NX, (i+1)*NX, NX, NX)   = Eigen::MatrixXd::Identity(NX, NX);
        Aeq.block(i*NX, N*NX + i*NU, NX, NU) = -B_d;
        beq.segment(i*NX, NX) = c_d;
    }

    // === Variable bounds (steering) ===
    Eigen::VectorXd lb = -1e5 * Eigen::VectorXd::Ones(nVar);
    Eigen::VectorXd ub =  1e5 * Eigen::VectorXd::Ones(nVar);
    for (int i = 0; i < N; ++i) {
        int idx = N * NX + i * NU;
        lb(idx) = -1.5; // min steering [rad]
        ub(idx) =  1.5; // max steering [rad]
    }

    // === Prepare qpOASES arrays ===
    real_t H_qp[nVar*nVar], g_qp[nVar], A_qp[nConstraints*nVar], lbA[nConstraints], ubA[nConstraints];

    // Copy Eigen matrices into qpOASES arrays using Eigen::Map
    Eigen::Map<Eigen::Matrix<real_t, nVar, nVar, Eigen::RowMajor>>(H_qp, nVar, nVar) = H;
    Eigen::Map<Eigen::Matrix<real_t, nVar, 1>>(g_qp, nVar) = g;
    Eigen::Map<Eigen::Matrix<real_t, nConstraints, nVar, Eigen::RowMajor>>(A_qp, nConstraints, nVar) = Aeq;
    Eigen::Map<Eigen::Matrix<real_t, nConstraints, 1>>(lbA, nConstraints) = beq;
    Eigen::Map<Eigen::Matrix<real_t, nConstraints, 1>>(ubA, nConstraints) = beq;
      

    // === Solve QP ===
    QProblem qp(nVar, nConstraints);
    Options options;
    options.setToMPC();
    qp.setOptions(options);

    int nWSR = 200;
    returnValue result = qp.init(H_qp, g_qp, A_qp, lb.data(), ub.data(), lbA, ubA, nWSR);
    
    //Initialize TargetState
    utfr_msgs::msg::TargetState target;

    if (result == SUCCESSFUL_RETURN) {
        real_t z[nVar];
        qp.getPrimalSolution(z);  
        double steering = z[N * NX]; // first control input (steering)
        steering = std::clamp(steering, -0.34, 0.34); //safety clamp max steering
        target.steering_angle = steering;
        setPrevD(steering);
        RCLCPP_INFO(this->get_logger(), "Optimal steering = %f", steering);
    } else {
          RCLCPP_WARN(this->get_logger(), "QP failed: %d", result);
          RCLCPP_ERROR(rclcpp::get_logger("MPC"), "QP solver failed with code %d", result);
    }

    splitVelocity(target, velocity_profile);

    return target;

}

void ControllerNode::splitVelocity(utfr_msgs::msg::TargetState &target, utfr_msgs::msg::VelocityProfile &velocity_profile){
    // Vehicle parameters config
    const double l_F = 0.811;  // a_cg - distance from CG to front axle [m]
    const double l_R = 0.719;  // b_cg - distance from CG to rear axle [m]
    const double w_f = 1.250 / 2.0;  // half front track width [m]
    const double w_r = 1.173 / 2.0;  // half rear track width [m]

    double delta = target.steering_angle;

    // Desired longitudinal velocity
    double v_desired = velocity_profile.velocities[20];

    // Bicycle model: calculate slip angle at CG
    double beta = atan((l_R / (l_F + l_R)) * tan(delta));
    double v_x_ = v_desired * cos(beta);
    double v_y_ = v_desired * sin(beta);
    double r_z = (v_x_ / (l_F + l_R)) * tan(delta);

    // Calculate wheel velocities
    auto wheelVelocity = [&](double x_offset, double y_offset) {
        double vx = v_x_ - r_z * y_offset;
        double vy = v_y_ + r_z * x_offset;
        return std::make_pair(vx, vy);
    };

    auto [vx_fl, vy_fl] = wheelVelocity(l_F, +w_f);   // front left
    auto [vx_fr, vy_fr] = wheelVelocity(l_F, -w_f);   // front right
    auto [vx_rl, vy_rl] = wheelVelocity(-l_R, +w_r);  // rear left
    auto [vx_rr, vy_rr] = wheelVelocity(-l_R, -w_r);  // rear right

    // Set target state values.
    target.speed_fl = vx_fl * cos(delta) + vy_fl * sin(delta);
    target.speed_fr = vx_fr * cos(delta) + vy_fr * sin(delta);
    target.speed_bl = vx_rl;
    target.speed_br = vx_rr;
} 
