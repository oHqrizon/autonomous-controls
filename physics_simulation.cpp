// Load Model

// ROS Publishers

// ROS Services

// ROS Subscriptions

// Init Params

// initVehicleModel

void RaceCarModelPlugin::initNoise(const sdf::ElementPtr &sdf) {
  std::string yaml_name = "";
  if (!sdf->HasElement("noise_config")) {
    RCLCPP_FATAL(_rosnode->get_logger(),
                 "gazebo_ros_race_car_model plugin missing <noise_config>, cannot proceed");
    return;
  } else {
    yaml_name = sdf->GetElement("noise_config")->Get<std::string>();
  }

  // Create noise object
  _noise = std::make_unique<eufs::models::Noise>(yaml_name);
}

void RaceCarModelPlugin::setPositionFromWorld() {
  _offset = _model->WorldPose();

  RCLCPP_DEBUG(_rosnode->get_logger(), "Got starting offset %f %f %f", _offset.Pos()[0],
               _offset.Pos()[1], _offset.Pos()[2]);

  if (!_random_spawn) {
    _state.x = 0.0;
    _state.y = 0.0;
    _state.z = 0.0;
    _state.yaw = 0.0;
  } else {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> randAngle(-M_PI / 8, M_PI / 8);
    std::uniform_real_distribution<double> randPos(-0.5, 0.5);
    _state.x = randPos(mt);
    _state.y = randPos(mt);
    _state.z = 0;
    _state.yaw = randAngle(mt);
  }
  _state.v_x = 0.0;
  _state.v_y = 0.0;
  _state.v_z = 0.0;
  _state.r_x = 0.0;
  _state.r_y = 0.0;
  _state.r_z = 0.0;
  _state.a_x = 0.0;
  _state.a_y = 0.0;
  _state.a_z = 0.0;
}

bool RaceCarModelPlugin::resetVehiclePosition(
    std::shared_ptr<std_srvs::srv::Trigger::Request>,
    std::shared_ptr<std_srvs::srv::Trigger::Response> response) {
  _state.x = 0.0;
  _state.y = 0.0;
  _state.z = 0.0;
  _state.yaw = 0.0;
  _state.v_x = 0.0;
  _state.v_y = 0.0;
  _state.v_z = 0.0;
  _state.r_x = 0.0;
  _state.r_y = 0.0;
  _state.r_z = 0.0;
  _state.a_x = 0.0;
  _state.a_y = 0.0;
  _state.a_z = 0.0;

  const ignition::math::Vector3d vel(0.0, 0.0, 0.0);
  const ignition::math::Vector3d angular(0.0, 0.0, 0.0);

  _model->SetWorldPose(_offset);
  _model->SetAngularVel(angular);
  _model->SetLinearVel(vel);

  return response->success;
}

void RaceCarModelPlugin::updateStateFromGazebo() {
    // Get world pose and velocities
    ignition::math::Pose3d pose = _model->WorldPose();
    ignition::math::Vector3d lin_vel_world = _model->WorldLinearVel();
    ignition::math::Vector3d ang_vel_world = _model->WorldAngularVel();

    // Compute offsets
    double dx = pose.Pos().X();
    double dy = pose.Pos().Y();
    double dz = pose.Pos().Z();
    double offset_yaw = _offset.Rot().Yaw();

    // Transform world position into planner/body frame
    _state.x =  dx * cos(-offset_yaw) - dy * sin(-offset_yaw);
    _state.y =  dx * sin(-offset_yaw) + dy * cos(-offset_yaw);
    _state.z =  dz;

    // Yaw relative to offset, wrapped
    _state.yaw = fmod(pose.Rot().Yaw() - offset_yaw + M_PI, 2.0 * M_PI) - M_PI;

    // Angular velocity in body frame (assuming small pitch/roll)
    _state.r_x = ang_vel_world.X();
    _state.r_y = ang_vel_world.Y();
    _state.r_z = ang_vel_world.Z();
}

// State to Car state msg

// Publish car state

// Publish wheel speeds

// Publish Odom + publish TF

// This runs at 1000hz
void RaceCarModelPlugin::updateState(const double dt) {
  if (!_command_Q.empty()) {
    gazebo::common::Time cmd_time = _cmd_time_Q.front();
    if ((_last_sim_time - cmd_time).Double() >= _control_delay) {
      std::shared_ptr<ackermann_msgs::msg::AckermannDriveStamped> cmd = _command_Q.front();
      _des_input.acc = cmd->drive.acceleration;
      _des_input.vel = cmd->drive.speed;
      _des_input.delta = cmd->drive.steering_angle;

      _command_Q.pop();
      _cmd_time_Q.pop();
    }
  }

  if (_command_mode == velocity) {
    double current_speed = std::sqrt(std::pow(_state.v_x, 2) + std::pow(_state.v_y, 2));
    _des_input.acc = (_des_input.vel - current_speed) / dt;
  }

  // If last command was more than 1s ago, then slow down car
  _act_input.acc = (_last_sim_time - _last_cmd_time) < 1.0 ? _des_input.acc : -1.0;
  // Make sure steering rate is within limits
  _act_input.delta +=
      (_des_input.delta - _act_input.delta >= 0 ? 1 : -1) *
      std::min(_max_steering_rate * dt, std::abs(_des_input.delta - _act_input.delta));
  _state.z = 0;

  _vehicle->updateState(_state, _act_input, dt);

  
  // ------------------ Parameters ------------------
  static double Kp_wheel = 500.0;       // wheel PI gain
  static double Ki_wheel = 50.0;
  static double Kv_wheel = 0.5;

  static double Kp_steer = 888.0;       // steering hinge PID (very nessecary 888 for goodluck)
  static double Kv_steer = 50.0;

  static double tau_max = 2000.0;       // max torque per wheel [Nm]
  static double tau_ramp_rate = 5000.0; // ramp-in
  static double tau_cap_runtime = 0.0;

  static double wheel_radius = _vehicle->getWheelRadius(); 

  // Integral terms for PI control (per wheel)
  static double ei_fl = 0, ei_fr = 0, ei_rl = 0, ei_rr = 0;

  // ------------------ Steering (front only) ------------------
  double delta_target = _act_input.delta;

  // Simple PD controller for hinge
  auto ctrl_steer = [&](gazebo::physics::JointPtr hinge){
      if (!hinge) return;
      double error = delta_target - hinge->Position(0);
      double torque = Kp_steer * error - Kv_steer * hinge->GetVelocity(0);
      hinge->SetForce(0, torque);
  };

  ctrl_steer(_left_steering_joint);
  ctrl_steer(_right_steering_joint);

  // ------------------ Wheel Targets ------------------
  //eufs_msgs::msg::WheelSpeeds wheel_speeds = _vehicle->getWheelSpeeds(_state, _act_input);

  // Convert RPM -> rad/s

  // ------------------ Wheel Targets ------------------
  double w_fl = 0.0, w_fr = 0.0, w_rl = 0.0, w_rr = 0.0;
  
  if (_state_machine->isManual()) {
      
      // --- Original Ackermann -> wheel speeds ---
      eufs_msgs::msg::WheelSpeeds wheel_speeds = _vehicle->getWheelSpeeds(_state, _act_input);

      // Convert RPM -> rad/s
      w_fl = wheel_speeds.lf_speed * 2.0 * M_PI / 60.0;
      w_fr = wheel_speeds.rf_speed * 2.0 * M_PI / 60.0;
      w_rl = wheel_speeds.lb_speed * 2.0 * M_PI / 60.0;
      w_rr = wheel_speeds.rb_speed * 2.0 * M_PI / 60.0;

  } else if (_last_target_state) {
      RCLCPP_INFO(rclcpp::get_logger("RaceCarModelPlugin"), "RUNNING CONTROLLER");
      // --- Controller TargetState -> wheel speeds ---
      w_fl = _last_target_state->speed_fl / wheel_radius;
      w_fr = _last_target_state->speed_fr / wheel_radius;
      w_rl = _last_target_state->speed_bl / wheel_radius;
      w_rr = _last_target_state->speed_br / wheel_radius;
  }

  // double w_fl = wheel_speeds.lf_speed * 2.0 * M_PI / 60.0;
  // double w_fr = wheel_speeds.rf_speed * 2.0 * M_PI / 60.0;
  // double w_rl = wheel_speeds.lb_speed * 2.0 * M_PI / 60.0;
  // double w_rr = wheel_speeds.rb_speed * 2.0 * M_PI / 60.0;

  // Ramp torque cap
  tau_cap_runtime = std::min(tau_max, tau_cap_runtime + tau_ramp_rate * dt);

  // Wheel PI controller
  auto ctrl_wheel = [&](gazebo::physics::JointPtr wheel_joint, double w_tgt, double &ei){
      if (!wheel_joint) return;
      double w_meas = wheel_joint->GetVelocity(0);
      double e = w_tgt - w_meas;
      ei += e * dt;
      double tau = Kp_wheel * e + Ki_wheel * ei - Kv_wheel * w_meas;
      tau = std::clamp(tau, -tau_cap_runtime, tau_cap_runtime);
      wheel_joint->SetForce(0, tau);
  };

  // ------------------ Apply Torque to Wheels ------------------
  // Only front or rear driven wheels get torque here
  ctrl_wheel(_left_front_wheel_joint,  w_fl, ei_fl);
  ctrl_wheel(_right_front_wheel_joint, w_fr, ei_fr);
  ctrl_wheel(_left_rear_wheel_joint,   w_rl, ei_rl); // rear wheels can be 0 if front-drive
  ctrl_wheel(_right_rear_wheel_joint,  w_rr, ei_rr);

  // ------------------ Update State from Gazebo ------------------ //
  updateStateFromGazebo();

  // ------------------ Publish ------------------
  double time_since_last_published = (_last_sim_time - _time_last_published).Double();
  if (time_since_last_published >= (1 / _publish_rate)) {
      _time_last_published = _last_sim_time;
      publishCarState();
      publishWheelSpeeds();
      publishOdom();
      publishTf();
  }

  // ------------------ State Machine ------------------
  _state_machine->spinOnce(_last_sim_time);
}


