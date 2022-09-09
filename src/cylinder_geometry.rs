use crate::{Float, PI};

// https://slideplayer.com/3188098/11/images/slide_1.jpg
pub struct CylinderGeometry {
  pub bore: Float,
  pub bore_area: Float,
  pub stroke: Float,
  pub compression: Float,
  pub displacement_volume: Float, // Vd = (UT - OT) x Area
  pub clearance_volume: Float, // Vc = OT
  pub volume: Float, // Vc + Vd
  pub con_rod_len: Float, // l
  pub crankshaft_radius: Float, // a. 2a = stroke
}

impl CylinderGeometry {
  pub fn from_data(bore: Float, stroke: Float, compression: Float) -> Self {
    let mut geom = CylinderGeometry {
      bore,
      stroke,
      compression,
      bore_area: 0.0,
      displacement_volume: 0.0,
      clearance_volume: 0.0,
      volume: 0.0,
      con_rod_len: 0.0,
      crankshaft_radius: stroke / 2.0,
    };

    geom.bore_area = PI * bore * bore / 4.0;
    geom.displacement_volume = geom.bore_area * stroke;
    geom.clearance_volume = geom.displacement_volume / (compression - 1.0);
    geom.volume = geom.displacement_volume + geom.clearance_volume;

    // 0.28 .. 0.33 https://de.wikipedia.org/wiki/Pleuelstangenverh%C3%A4ltnis
    geom.con_rod_len = geom.crankshaft_radius / 0.3;

    return geom;
  }

  pub fn piston_bolt_position_over_angle(&self, angle: Float) -> Float { // called "s"
    // a * cos angle + (l^2 - a^2 sin^2 angle) ^ 1/2
    // a is crankshaft diameter, l is pleuellänge

    let a = self.crankshaft_radius;
    let l = self.con_rod_len;

    let sine_angle = Float::sin(angle);
    let term = l * l - (a * a * sine_angle * sine_angle);
    return a * Float::cos(angle) + Float::powf(term, 1.0/2.0);
  }

  pub fn volume_angle(&self, piston_bolt_pos: Float) -> Float {
    self.clearance_volume + self.bore_area * (self.con_rod_len + self.crankshaft_radius - piston_bolt_pos)
  }
}