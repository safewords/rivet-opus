//! Opus (RFC 6716, as updated by RFC 8251), both ways, in Rust: no C, no
//! system libraries, no build script.
//!
//! - [`Decoder`]: one Opus stream (mono or stereo) to interleaved `f32` at
//!   8, 12, 16, 24 or 48 kHz — SILK, CELT and hybrid frames, every frame
//!   size and packet code, mode transitions with and without redundancy,
//!   packet loss concealment and LBRR forward error correction.
//! - [`Encoder`]: interleaved `f32` at those rates to packets — CELT,
//!   SILK and hybrid modes, 2.5 to 60 ms, 6 to 510 kb/s, CBR or VBR,
//!   optional LBRR.
//! - [`MultistreamDecoder`] / [`MultistreamEncoder`] and [`OpusHead`]:
//!   RFC 7845 channel mapping families 0 and 1 (mono to 7.1; family 255 on
//!   decode).
//! - [`packet`]: the RFC 6716 §3 framing, parsed and built.
//!
//! PCM is interleaved `f32`, ±1.0 full scale.
//!
//! ```
//! use opus::{Decoder, Encoder, EncoderConfig};
//!
//! let cfg = EncoderConfig { channels: 2, bitrate: 96_000, ..EncoderConfig::default() };
//! let mut enc = Encoder::new(cfg)?;
//! let mut dec = Decoder::new(48_000, 2)?;
//! let pcm = vec![0.0f32; enc.frame_samples() * 2];
//! let packet = enc.encode(&pcm)?;
//! let out = dec.decode(Some(&packet))?;
//! assert_eq!(out.len(), pcm.len());
//! // The decoded audio lags the input by enc.lookahead() samples at
//! // 48 kHz: the stream's pre-skip.
//! # Ok::<(), opus::Error>(())
//! ```

// Unsafe code is confined to the vector kernels of `simd`.
#![deny(unsafe_code)]
#![warn(missing_docs)]
// Signal processing indexes several arrays in step; the index loops read
// closer to the RFC's formulas than iterator chains would.
// Shift-and-add expressions are written as the RFC writes them.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments, clippy::precedence)]

mod celt;
mod decoder;
mod encoder;
mod error;
mod mdct;
mod multistream;
pub mod packet;
mod range;
mod resample;
mod silk;
#[allow(unsafe_code)]
mod simd;

pub use decoder::{Decoder, SAMPLE_RATES};
pub use encoder::{Application, Encoder, EncoderConfig, LOOKAHEAD_48K};
pub use error::{Error, Result};
pub use multistream::{MultistreamDecoder, MultistreamEncoder, OpusHead, family1_layout};
pub use packet::{Bandwidth, Mode};
pub use range::{RangeDecoder, RangeEncoder};
