mod cylinder_geometry;
mod engine_state;

type Float = f32;
const PI: Float = std::f32::consts::PI;
// type Float = f64;
// const PI: Float = std::f64::consts::PI;

const KAPPA_AIR: Float = 1.40;
const CV_AIR: Float = 0.718; // kJ / kg * K
const ATMOSPHERIC_PRESSURE: Float = 100000.0; // atmospheric, 1 bar, 100 kPa
const ATMOSPHERIC_TEMP_C: Float = 20.0;
const DENSITY_AFR: Float = 1.29; // kg / m^3

const CSV: bool = true;

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Write;
    use hound::SampleFormat;
    use crate::cylinder_geometry::CylinderGeometry;
    use crate::{Float, KAPPA_AIR, PI, ATMOSPHERIC_PRESSURE, ATMOSPHERIC_TEMP_C, DENSITY_AFR, CV_AIR, CSV};
    use crate::engine_state::MechanicalEngineState;

    #[test]
    fn cyl() {
        // https://www.carfolio.com/alfa-romeo-gtv6-2.5-28336
        let cyl = CylinderGeometry::from_data(88.0 / 1000.0, 68.3 / 1000.0, 9.0);

        // compression = (Vc + Vd) / Vc = 1 + Vd / Vc
        // compression - 1 = Vd/Vc
        // (compression - 1) / Vd = 1 / Vc
        // Vd / (compression - 1) = Vc

        // TODO: we can't allow a total vacuum at load = 0, right?
        let load = 1.0;
        let pressure_ut = load * ATMOSPHERIC_PRESSURE; // being off-load constructs a vacuum (very simplified formula)
        //let pressure_ot = pressure_ut * cyl.volume / cyl.clearance_volume; // boyle's law
        // We can't use Boyle's Law here because the air will also heat due to compression
        let pressure_ot = isentrope_pressure(pressure_ut, cyl.volume, cyl.clearance_volume, KAPPA_AIR);

        let cyl_angle_offset: Vec<Float> = vec![0.0, 360.0].into_iter().map(| x| x * PI / 180.0).collect();

        let mut csv = File::create("csv.csv").unwrap();
        write!(csv, "degree,cyl,pressure,ign,volume\n").unwrap();

        let spec = hound::WavSpec {
            channels: 1,
            sample_rate: 23760 * 2,//44100,
            bits_per_sample: 32,
            sample_format: SampleFormat::Float
        };

        let mut writer = hound::WavWriter::create("./2k.wav", spec).unwrap();

        // get the ign pressure at OT, because that's relevant for the expansion pressure level
        let burn_curve_ot = ignition_get_burn_curve(0.0);
        let ign_ot = pressure_by_ignition(DENSITY_AFR * cyl.volume, CV_AIR, 0.0, pressure_ut, burn_curve_ot);

        let mut samples = Vec::<f32>::new();
        for x in 0..720/* * 33 * 4 /* 33 rps = 2k rpm */ */ {
            let angle_radians: Float = x as Float * PI / 180.0;
            let mut pressure_sample = 0.0;

            for idx in 0..1 {
                //let cylinder_angle = angle_radians + cylinder_offset * (idx as Float);
                let cylinder_angle = (angle_radians + *cyl_angle_offset.get(idx).unwrap()) % (4.0 * PI);
                let state = MechanicalEngineState::from_radians(cylinder_angle);
                let bolt_pos = cyl.piston_bolt_position_over_angle(cylinder_angle);
                let volume = cyl.volume_angle(bolt_pos);

                if state == MechanicalEngineState::COMPRESSION {
                    //let pressure = pressure_ut * volume / cyl.clearance_volume;
                    let mut pressure = isentrope_pressure(pressure_ut, cyl.volume, volume, KAPPA_AIR);
                    // For now, let's assume a linear burning.
                    let burn_curve = ignition_get_burn_curve(cylinder_angle - 2.0 * PI);
                    let mass = DENSITY_AFR * volume;
                    let pi = pressure_by_ignition(mass, CV_AIR, 0.0, pressure_ut, burn_curve);
                    pressure += pi;
                    //pressure += fake_pressure_by_ignition(cylinder_angle - 2.0 * PI) * ATMOSPHERIC_PRESSURE * 10.0 * load;
                    //println!("[{idx}]: Volume at {x} degrees ({bolt_pos} m): {volume}. Pressure: {pressure}");
                    if CSV {
                        write!(csv, "{x},{},{pressure},{pi},{volume}\n", idx + 1).unwrap();
                    }
                    pressure_sample += pressure;
                } else if state == MechanicalEngineState::WORK {
                    // For now, let's assume a linear burning.
                    let mass = DENSITY_AFR * volume;
                    let burn_curve = ignition_get_burn_curve(cylinder_angle - 2.0 * PI);

                    // expansion level reduced BY ignition pressure at OT, because the actual pressure is added later on
                    let mut pressure = isentrope_pressure(pressure_ot /* + ign_ot*/, cyl.clearance_volume, volume, KAPPA_AIR);// - ign_ot;
                    let ign = pressure_by_ignition(mass, CV_AIR, 0.0, pressure_ut, burn_curve);
                    pressure += ign;
                    //pressure += fake_pressure_by_ignition(cylinder_angle - 2.0 * PI) * ATMOSPHERIC_PRESSURE * 10.0 * load;

                    if CSV {
                        write!(csv, "{x},{},{pressure},{ign},{volume}\n", idx + 1).unwrap();
                    }
                    pressure_sample += pressure;
                } else if state == MechanicalEngineState::EXHAUST {
                    // How about pressure_ut + burn energy?
                    let mass = DENSITY_AFR * volume;
                    let mut ign = pressure_by_ignition(mass, CV_AIR, 0.0, pressure_ut, 1.0);

                    // maybe let ign deplete over the whole angle? [3Pi, 4Pi]
                    let amnt : Float = (cylinder_angle / PI) - 3.0;
                    ign *= 1.0 - (amnt / 2.0);
                    let pressure = pressure_ut + ign;

                    if CSV {
                        write!(csv, "{x},{},{pressure},{ign},{volume}\n", idx + 1).unwrap();
                    }
                    pressure_sample += pressure;
                } else { // INTAKE
                    // when intaking, pressure starts at atmospheric and volume is increased, so pressure gently sinks.
                    // let pressure = 100000.0 * cyl.clearance_volume / volume;

                    // The gleichraumprozess actually does both exhausting and intaking together.
                    //let pressure = isochoric_pressure(pressure_ot, 275.0 + 1200.0, 275.0 + 20.0);

                    // 3rd try: the pressure stays roughly the same, but the pressure niveau depends on the load. A load of 1.0 is just below atmospheric pressure.
                    let pressure = pressure_ut;
                    if CSV {
                        write!(csv, "{x},{},{pressure},{0},{volume}\n", idx + 1).unwrap();
                    }
                    pressure_sample += pressure;

                }
            }
            //pressure_sample /= 6.0;
            if CSV {
                write!(csv, "{x},0,{pressure_sample},0,0\n").unwrap();
            }
            //writer.write_sample(((pressure_sample / pressure_ot) - 0.5) as f32).unwrap();
            //writer.write_sample(((pressure_sample / (6.0 * (pressure_ot + ign_ot))) - 0.5) as f32).unwrap();
            //writer.write_sample(pressure_sample).unwrap();
            samples.push(pressure_sample);
        }

        let max = *samples.iter().max_by(|x, y| x.abs().partial_cmp(&y.abs()).unwrap()).unwrap();
        for sample in samples {
            writer.write_sample((sample / max) - 0.5).unwrap();
        }

        writer.finalize().unwrap();
    }

    // the problem with the "Gleichraumprozess" is, it isn't real and matching the 4 strokes?

    /// Kappa is also known as "isentropen exponent" (heat capacity ratio)
    #[inline]
    fn isentrope_pressure(p0: Float, v0: Float, v: Float, kappa: Float) -> Float {
        p0 * ((v0 / v).powf(kappa))
    }

    #[inline]
    fn isentrope_temperature(t0: Float, v0: Float, v: Float, kappa: Float) -> Float {
        t0 * ((v0 / v).powf(kappa - 1.0))
    }

    // p/T == const
    #[inline]
    fn isochoric_pressure(p0: Float, t0: Float, t1: Float) -> Float {
        t1/t0 * p0
    }

    // The energy/enthalpy released by ignition is increasing the internal energy of the system.
    // This temperature increase will then be distribute like an isochoric temperature increase.

    #[inline]
    fn fake_pressure_by_ignition(angle_relative_to_ot: Float) -> Float {
        const ADVANCE: Float = 10.0;
        const ADVANCE_RADIANS: Float = ADVANCE * PI / 180.0;

        if angle_relative_to_ot.abs() > ADVANCE_RADIANS {
            0.0
        } else {
            let x = angle_relative_to_ot / ADVANCE_RADIANS;
            //-(x * x) + 1.0 // parabola where { -1, 1 } are zero crossings and the max is at (0, 1)
            -1.0/2.0 * ((x - 0.5) * (x - 0.5)) + 1.0
        }
    }

    #[inline]
    fn energy_by_ignition(mass: Float, _afr: Float) -> Float {
        44000.0 * mass / 15.0 // 44k kJ = 44 Mj. / 15 as a rough way to come back to fuel mass
    }

    /// deltaQ = m * c_v * deltaT. Since c_v is kJ / kg * K, deltaQ needs to be in kJ
    #[inline]
    fn temperature_rise_by_energy(delta_q: Float, mass: Float, c_v: Float) -> Float {
        delta_q / (mass * c_v)
    }

    // TODO: Should p0 be ATMOSPHERIC or rather the compressed pressure up to OT? Something?
    #[inline]
    fn pressure_by_ignition(mass: Float, c_v: Float, afr: Float, p0: Float, burn_progress: Float) -> Float {
        if burn_progress == 0.0 {
            return 0.0;
        }

        let delta_q = energy_by_ignition(mass, afr) * burn_progress;
        let t0 = 273.15 + ATMOSPHERIC_TEMP_C;
        let t1 = t0 + temperature_rise_by_energy(delta_q, mass, c_v);
        isochoric_pressure(p0, t0, t1)
    }

    /// Linear for now
    #[inline]
    fn ignition_get_burn_curve(angle_relative_to_ot: Float) -> Float {
        // const ADVANCE: Float = 40.0;
        // const ADVANCE_RADIANS: Float = ADVANCE * PI / 180.0;
        // if angle_relative_to_ot < -ADVANCE_RADIANS {
        //     0.0
        // } else if angle_relative_to_ot > ADVANCE_RADIANS {
        //     1.0
        // } else {
        //     let t = angle_relative_to_ot / ADVANCE_RADIANS;
        //     if t < 0.0 {
        //         0.0
        //     } else {
        //         t
        //     }
        // }

        const ADVANCE: Float = 20.0;
        const WIDTH: Float = 60.0;
        const ADVANCE_RADIANS: Float = ADVANCE * PI / 180.0;
        const WIDTH_RADIANS: Float = WIDTH * PI / 180.0;

        if angle_relative_to_ot < -ADVANCE_RADIANS {
            0.0
        } else if angle_relative_to_ot > WIDTH_RADIANS /*- ADVANCE_RADIANS */{
            1.0
        } else {
            let x = (angle_relative_to_ot + ADVANCE_RADIANS) / (WIDTH_RADIANS + ADVANCE_RADIANS);
            x
        }
    }
}
