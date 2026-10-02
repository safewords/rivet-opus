//! Not a test: a tool, run by hand (`cargo test --test fit_resampler --
//! --ignored --nocapture` with OPUS_TESTVECTORS set), that identifies the
//! response of the reference decoder's (non-normative) SILK resampler from
//! the official test vectors, treated as black-box data: for SILK-only
//! vectors at one bandwidth it decodes the SILK signal at its internal rate
//! with this crate, and fits, by least squares, the causal polyphase filter
//! that best maps it to the reference's 48 kHz output.
use std::path::PathBuf;

fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for i in 0..n {
        let piv = (i..n).max_by(|&x, &y| a[x][i].abs().total_cmp(&a[y][i].abs())).unwrap();
        a.swap(i, piv);
        b.swap(i, piv);
        for r in i + 1..n {
            let f = a[r][i] / a[i][i];
            for c in i..n {
                a[r][c] -= f * a[i][c];
            }
            b[r] -= f * b[i];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for c in i + 1..n {
            s -= a[i][c] * x[c];
        }
        x[i] = s / a[i][i];
    }
    x
}

#[test]
#[ignore]
fn fit() {
    let d = PathBuf::from(std::env::var("OPUS_TESTVECTORS").unwrap());
    let k: usize = std::env::var("TAPS").ok().and_then(|v| v.parse().ok()).unwrap_or(24);
    for (v, rate) in [(2usize, 8000u32), (3, 12000), (4, 16000)] {
        let l = (48000 / rate) as usize;
        let b = std::fs::read(d.join(format!("testvector{v:02}.bit"))).unwrap();
        let r: Vec<f64> = std::fs::read(d.join(format!("testvector{v:02}.dec"))).unwrap().chunks_exact(2).map(|c| f64::from(i16::from_le_bytes([c[0], c[1]])) / 32768.0).collect();
        let mut dec = opus::Decoder::new(rate, 2).unwrap();
        let mut x: Vec<f64> = Vec::new();
        let mut p = 0;
        while p + 8 <= b.len() {
            let len = u32::from_be_bytes(b[p..p + 4].try_into().unwrap()) as usize;
            p += 8;
            x.extend(dec.decode(Some(&b[p..p + len])).unwrap().iter().map(|&s| f64::from(s)));
            p += len;
        }
        let m_count = (x.len() / 2).min(r.len() / 2 / l);
        let mut total_s = 0.0;
        let mut total_e = 0.0;
        let mut out = format!("    // {rate} Hz -> 48 kHz, {k} taps per phase, phase-major, newest input first\n    &[\n");
        for ph in 0..l {
            let mut ata = vec![vec![0.0; k]; k];
            let mut aty = vec![0.0; k];
            for m in k..m_count {
                for c in 0..2 {
                    let y = r[2 * (l * m + ph) + c];
                    for i in 0..k {
                        let xi = x[2 * (m - i) + c];
                        aty[i] += xi * y;
                        for j in 0..k {
                            ata[i][j] += xi * x[2 * (m - j) + c];
                        }
                    }
                }
            }
            for i in 0..k {
                ata[i][i] += 1e-9;
            }
            let g = solve(ata, aty);
            for m in k..m_count {
                for c in 0..2 {
                    let y = r[2 * (l * m + ph) + c];
                    let est: f64 = (0..k).map(|i| g[i] * x[2 * (m - i) + c]).sum();
                    total_s += y * y;
                    total_e += (y - est) * (y - est);
                }
            }
            out.push_str("        &[");
            out.push_str(&g.iter().map(|v| format!("{:.9}", v)).collect::<Vec<_>>().join(", "));
            out.push_str("],\n");
        }
        out.push_str("    ],\n");
        eprintln!("vector {v} ({rate} Hz): fit SNR {:.2} dB", 10.0 * (total_s / total_e).log10());
        if std::env::var_os("PRINT").is_some() {
            eprintln!("{out}");
        }
    }
}
