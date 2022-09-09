use crate::{Float, PI};

#[derive(PartialEq)]
pub enum MechanicalEngineState {
  INTAKE,
  COMPRESSION,
  WORK,
  EXHAUST
}

impl MechanicalEngineState {
  pub fn from_radians(angle: Float) -> MechanicalEngineState {
    let modulo = angle % (4.0 * PI);
    let index = modulo / PI;

    if index <= 1.0 {
      return MechanicalEngineState::INTAKE;
    } else if index <= 2.0 {
      return MechanicalEngineState::COMPRESSION;
    } else if index <= 3.0 {
      return MechanicalEngineState::WORK;
    } else if index <= 4.0 {
      return MechanicalEngineState::EXHAUST;
    } else {
      panic!("Unknown Engine State: {}", index);
    }
  }
}