//! Opus (RFC 6716, as updated by RFC 8251), both ways, in Rust.

#![forbid(unsafe_code)]

mod celt;
mod decoder;
mod error;
mod mdct;
pub mod packet;
mod range;
mod resample;
mod silk;

pub use decoder::{Decoder, SAMPLE_RATES};
pub use error::{Error, Result};
pub use packet::{Bandwidth, Mode};
pub use range::{RangeDecoder, RangeEncoder};
