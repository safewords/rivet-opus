//! Opus (RFC 6716, as updated by RFC 8251), both ways, in Rust.

#![forbid(unsafe_code)]

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
