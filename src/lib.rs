//! Opus (RFC 6716, as updated by RFC 8251), both ways, in Rust.

#![forbid(unsafe_code)]

mod celt;
mod error;
mod mdct;
pub mod packet;
mod range;

pub use error::{Error, Result};
pub use range::{RangeDecoder, RangeEncoder};
