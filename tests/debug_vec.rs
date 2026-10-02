use std::path::PathBuf;
#[test]
#[ignore]
fn debug_mismatch() {
    let d = PathBuf::from(std::env::var("OPUS_TESTVECTORS").unwrap());
    let v: usize = std::env::var("VEC").unwrap().parse().unwrap();
    let b = std::fs::read(d.join(format!("testvector{v:02}.bit"))).unwrap();
    let mut p = 0;
    let mut dec = opus::Decoder::new(48000, 2).unwrap();
    let mut i = 0;
    let mut shown = 0;
    let mut stats = std::collections::BTreeMap::new();
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        let range = u32::from_be_bytes(b[p + 4..p + 8].try_into().unwrap());
        p += 8;
        let data = &b[p..p + len];
        p += len;
        dec.decode(Some(data)).unwrap();
        let toc = data[0];
        let key = (toc >> 3, (toc >> 2) & 1, toc & 3);
        let e = stats.entry(key).or_insert((0, 0));
        e.0 += 1;
        eprintln!("== packet {i} ok {}", dec.final_range() == range);
        if dec.final_range() != range {
            e.1 += 1;
            if shown < 10 {
                eprintln!("packet {i}: config {} stereo {} code {} len {len}", toc >> 3, (toc >> 2) & 1, toc & 3);
                shown += 1;
            }
        }
        i += 1;
    }
    for (k, v) in stats {
        eprintln!("{k:?}: {} packets, {} bad", v.0, v.1);
    }
}

#[test]
#[ignore]
fn debug_align() {
    let d = PathBuf::from(std::env::var("OPUS_TESTVECTORS").unwrap());
    let v: usize = std::env::var("VEC").unwrap().parse().unwrap();
    let b = std::fs::read(d.join(format!("testvector{v:02}.bit"))).unwrap();
    let r: Vec<f64> = std::fs::read(d.join(format!("testvector{v:02}.dec"))).unwrap().chunks_exact(2).map(|c| f64::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0).collect();
    let mut p = 0;
    let mut dec = opus::Decoder::new(48000, 2).unwrap();
    let mut out: Vec<f64> = Vec::new();
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        p += 8;
        out.extend(dec.decode(Some(&b[p..p + len])).unwrap().iter().map(|&x| f64::from(x)));
        p += len;
    }
    eprintln!("lengths ours {} ref {}", out.len(), r.len());
    // Best lag in samples (per channel) and gain over a middle stretch.
    let n = out.len().min(r.len()) / 2;
    let (a, z) = (n / 4, n / 4 + 48000 * 4);
    for lag in -40i64..=40 {
        let (mut xy, mut yy, mut xx) = (0.0, 0.0, 0.0);
        for i in a..z {
            let j = (i as i64 + lag) as usize;
            for c in 0..2 {
                xy += r[2 * i + c] * out[2 * j + c];
                yy += out[2 * j + c] * out[2 * j + c];
                xx += r[2 * i + c] * r[2 * i + c];
            }
        }
        let g = xy / yy;
        let err = xx - 2.0 * g * xy + g * g * yy;
        eprintln!("lag {lag}: gain {g:.4} snr {:.2}", 10.0 * (xx / err).log10());
    }
}

#[test]
#[ignore]
fn debug_lowpass() {
    let d = PathBuf::from(std::env::var("OPUS_TESTVECTORS").unwrap());
    let v: usize = std::env::var("VEC").unwrap().parse().unwrap();
    let b = std::fs::read(d.join(format!("testvector{v:02}.bit"))).unwrap();
    let r: Vec<f64> = std::fs::read(d.join(format!("testvector{v:02}.dec"))).unwrap().chunks_exact(2).map(|c| f64::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0).collect();
    let mut p = 0;
    let mut dec = opus::Decoder::new(48000, 2).unwrap();
    let mut out: Vec<f64> = Vec::new();
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        p += 8;
        out.extend(dec.decode(Some(&b[p..p + len])).unwrap().iter().map(|&x| f64::from(x)));
        p += len;
    }
    for w in [1usize, 4, 8, 16, 32, 64] {
        // Moving average of w samples on the left channel.
        let ma = |x: &[f64]| -> Vec<f64> {
            let n = x.len() / 2;
            let mut acc = 0.0;
            let mut o = vec![0.0; n];
            for i in 0..n {
                acc += x[2 * i];
                if i >= w { acc -= x[2 * (i - w)]; }
                o[i] = acc / w as f64;
            }
            o
        };
        let (a, bb) = (ma(&r), ma(&out));
        let (mut s, mut e) = (0.0, 0.0);
        for i in 0..a.len().min(bb.len()) { s += a[i] * a[i]; e += (a[i] - bb[i]).powi(2); }
        eprintln!("ma {w}: snr {:.2}", 10.0 * (s / e).log10());
    }
}

#[test]
#[ignore]
fn debug_transfer() {
    let d = PathBuf::from(std::env::var("OPUS_TESTVECTORS").unwrap());
    let v: usize = std::env::var("VEC").unwrap().parse().unwrap();
    let b = std::fs::read(d.join(format!("testvector{v:02}.bit"))).unwrap();
    let r: Vec<f64> = std::fs::read(d.join(format!("testvector{v:02}.dec"))).unwrap().chunks_exact(2).map(|c| f64::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0).collect();
    let mut p = 0;
    let mut dec = opus::Decoder::new(48000, 2).unwrap();
    let mut out: Vec<f64> = Vec::new();
    while p + 8 <= b.len() {
        let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
        p += 8;
        out.extend(dec.decode(Some(&b[p..p + len])).unwrap().iter().map(|&x| f64::from(x)));
        p += len;
    }
    let n = out.len().min(r.len()) / 2;
    let fl = 2048;
    for f in [200.0, 500.0, 1000.0, 2000.0, 3000.0, 3500.0, 4500.0, 6000.0, 7000.0, 7500.0, 9000.0, 11000.0, 14000.0, 18000.0] {
        let w = 2.0 * std::f64::consts::PI * f / 48000.0;
        let (mut cre, mut cim, mut po, mut pr) = (0.0, 0.0, 0.0, 0.0);
        let mut s = 0;
        while s + fl < n {
            let (mut rr, mut ri, mut or, mut oi) = (0.0, 0.0, 0.0, 0.0);
            for k in 0..fl {
                let win = 0.5 - 0.5 * (2.0 * std::f64::consts::PI * k as f64 / fl as f64).cos();
                let (c, sn) = ((w * k as f64).cos(), (w * k as f64).sin());
                let x = r[2 * (s + k)] * win;
                let y = out[2 * (s + k)] * win;
                rr += x * c; ri -= x * sn; or += y * c; oi -= y * sn;
            }
            // R * conj(O)
            cre += rr * or + ri * oi;
            cim += ri * or - rr * oi;
            po += or * or + oi * oi;
            pr += rr * rr + ri * ri;
            s += fl;
        }
        let mag = (cre * cre + cim * cim).sqrt() / po;
        let ph = cim.atan2(cre);
        let coh = (cre * cre + cim * cim) / (po * pr);
        eprintln!("f {f}: |H| {mag:.3} phase {ph:.3} delay {:.3} samples coherence {coh:.4} level {:.1} dB", -ph / w, 10.0 * (pr / 1e3).log10());
    }
}
