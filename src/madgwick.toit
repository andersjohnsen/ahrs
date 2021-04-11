import math
import math3d

class Madgwick:
  static BETA_DEFAULT_/float ::= 0.1

  beta_/float := BETA_DEFAULT_
  rotation_/math3d.Quaternion := math3d.Quaternion.identity

  constructor:

  rotation -> math3d.Quaternion: return rotation_

  update_imu gyro/math3d.Vector3 accel/math3d.Vector3 elapsed/Duration:
    elapsed_s := elapsed.in_ns.to_float / Duration.NANOSECONDS_PER_SECOND

    r := rotation_

    q1 := r.w
    q2 := r.x
    q3 := r.y
    q4 := r.z

    gx := gyro.x
    gy := gyro.y
    gz := gyro.z

    _2q1 := 2.0 * q1
    _2q2 := 2.0 * q2
    _2q3 := 2.0 * q3
    _2q4 := 2.0 * q4
    _4q1 := 4.0 * q1
    _4q2 := 4.0 * q2
    _4q3 := 4.0 * q3
    _8q2 := 8.0 * q2
    _8q3 := 8.0 * q3
    q1q1 := q1 * q1
    q2q2 := q2 * q2
    q3q3 := q3 * q3
    q4q4 := q4 * q4

    a_length := accel.length

    // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
    if a_length == 0.0: return

    // Normalise accelerometer measurement
    ax := accel.x / a_length
    ay := accel.y / a_length
    az := accel.z / a_length

    s1 := _4q1 * q3q3 + _2q3 * ax + _4q1 * q2q2 - _2q2 * ay
    s2 := _4q2 * q4q4 - _2q4 * ax + 4.0 * q1q1 * q2 - _2q1 * ay - _4q2 + _8q2 * q2q2 + _8q2 * q3q3 + _4q2 * az
    s3 := 4.0 * q1q1 * q3 + _2q1 * ax + _4q3 * q4q4 - _2q4 * ay - _4q3 + _8q3 * q2q2 + _8q3 * q3q3 + _4q3 * az
    s4 := 4.0 * q2q2 * q4 - _2q2 * ax + 4.0 * q3q3 * q4 - _2q3 * ay
    norm := 1.0 / (math.sqrt s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4)    // normalise step magnitude
    s1 *= norm;
    s2 *= norm;
    s3 *= norm;
    s4 *= norm;

    // Compute rate of change of quaternion.
    q_dot1 := 0.5 * (-q2 * gx - q3 * gy - q4 * gz) - beta_ * s1
    q_dot2 := 0.5 * (q1 * gx + q3 * gz - q4 * gy) - beta_ * s2
    q_dot3 := 0.5 * (q1 * gy - q2 * gz + q4 * gx) - beta_ * s3
    q_dot4 := 0.5 * (q1 * gz + q2 * gy - q3 * gx) - beta_ * s4

    q1 += q_dot1 * elapsed_s
    q2 += q_dot2 * elapsed_s
    q3 += q_dot3 * elapsed_s
    q4 += q_dot4 * elapsed_s

    norm = 1.0 / (math.sqrt q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4)

    rotation_ = math3d.Quaternion
      q2 * norm
      q3 * norm
      q4 * norm
      q1 * norm
