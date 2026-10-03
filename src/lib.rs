//! Opus (RFC 6716, as updated by RFC 8251), both ways, in Rust.

#![forbid(unsafe_code)]
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
mod resample_fit;
mod silk;

pub use decoder::{Decoder, SAMPLE_RATES};
pub use encoder::{Application, Encoder, EncoderConfig, LOOKAHEAD_48K};
pub use error::{Error, Result};
pub use multistream::{MultistreamDecoder, MultistreamEncoder, OpusHead, family1_layout};
pub use packet::{Bandwidth, Mode};
pub use range::{RangeDecoder, RangeEncoder};
