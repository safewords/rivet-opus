//! The crate's one error type.

/// What went wrong. Every malformed input comes back as one of these; the
/// decoder never panics on bytes it is given.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The packet breaks the framing rules of RFC 6716 §3 (\[R1\]–\[R7\]) or a
    /// frame's contents are not decodable.
    #[error("invalid Opus packet: {0}")]
    InvalidPacket(String),
    /// An `OpusHead` (RFC 7845 §5.1) that cannot be used.
    #[error("invalid OpusHead: {0}")]
    InvalidHead(String),
    /// Valid Opus this crate does not implement, named.
    #[error("unsupported Opus feature: {0}")]
    Unsupported(String),
    /// A configuration the caller asked for that cannot be used: a sample
    /// rate, channel count, bit rate or frame size outside what Opus allows.
    #[error("invalid Opus configuration: {0}")]
    Config(String),
    /// The caller's buffer or input has the wrong size.
    #[error("bad argument: {0}")]
    BadArgument(String),
}

pub(crate) fn invalid(msg: impl Into<String>) -> Error {
    Error::InvalidPacket(msg.into())
}

pub(crate) fn config(msg: impl Into<String>) -> Error {
    Error::Config(msg.into())
}

/// The crate's `Result`.
pub type Result<T> = std::result::Result<T, Error>;
