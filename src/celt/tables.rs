//! Constants and data tables of the CELT layer (the 48 kHz mode of RFC 6716
//! §4.3). The arrays below are data copied from the output of
//! `tools/appendix_a_tables.py` (RFC 6716 Appendix A read as data); the
//! tables that CELT_SPEC derives (pulse cache, caps, `LOG2_FRAC_TABLE`,
//! `logN`, interleaving orders, `exp2_table8`) are computed in
//! [`super::mode`] instead.

/// The overlap of the low-overlap window: 2.5 ms at 48 kHz (CELT_SPEC §10.4).
pub const OVERLAP: usize = 120;

/// Number of energy bands (CELT_SPEC conventions; RFC 6716 Table 55).
pub const NB_EBANDS: usize = 21;

/// Allocation quantities are in 1/2^BITRES bit (CELT_SPEC conventions).
pub const BITRES: i32 = 3;

/// Most fine-energy bits per band and channel (CELT_SPEC §3.2, §13).
pub const MAX_FINE_BITS: i32 = 8;

/// The fine-energy offset of the allocation split (CELT_SPEC §6.9).
pub const FINE_OFFSET: i32 = 21;

/// Theta resolution offsets (CELT_SPEC §8.2.4).
pub const QTHETA_OFFSET: i32 = 4;
pub const QTHETA_OFFSET_TWOPHASE: i32 = 16;

/// Largest pseudo-pulse index of the pulse cache (CELT_SPEC §7.2).
pub const MAX_PSEUDO: usize = 40;

/// Band edges at LM = 0 (`EBAND5MS` as `usize`).
pub const EBANDS: [usize; NB_EBANDS + 1] = {
    let mut e = [0usize; NB_EBANDS + 1];
    let mut i = 0;
    while i <= NB_EBANDS {
        e[i] = EBAND5MS[i] as usize;
        i += 1;
    }
    e
};

/// The pre-emphasis / de-emphasis coefficient (RFC 6716 §4.3.7.2,
/// CELT_SPEC §10.6; `PREEMPH[0]`).
pub const PREEMPH_COEF: f32 = PREEMPH[0];

/// RFC 6716 Appendix A celt/static_modes_float.h, `mode48000_960_120.preemph`, extracted by tools/appendix_a_tables.py.
pub const PREEMPH: [f32; 4] = [0.8500061, 0.0, 1.0, 1.0];

/// RFC 6716 Appendix A celt/modes.c, `eband5ms`, extracted by tools/appendix_a_tables.py.
pub const EBAND5MS: [i16; 22] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 34, 40, 48, 60, 78, 100,
];

/// RFC 6716 Appendix A celt/modes.c, `band_allocation`, extracted by tools/appendix_a_tables.py.
pub const BAND_ALLOCATION: [[u8; 21]; 11] = [
    [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ],
    [
        90, 80, 75, 69, 63, 56, 49, 40, 34, 29, 20, 18, 10, 0, 0, 0, 0, 0, 0, 0, 0,
    ],
    [
        110, 100, 90, 84, 78, 71, 65, 58, 51, 45, 39, 32, 26, 20, 12, 0, 0, 0, 0, 0, 0,
    ],
    [
        118, 110, 103, 93, 86, 80, 75, 70, 65, 59, 53, 47, 40, 31, 23, 15, 4, 0, 0, 0, 0,
    ],
    [
        126, 119, 112, 104, 95, 89, 83, 78, 72, 66, 60, 54, 47, 39, 32, 25, 17, 12, 1, 0, 0,
    ],
    [
        134, 127, 120, 114, 103, 97, 91, 85, 78, 72, 66, 60, 54, 47, 41, 35, 29, 23, 16, 10, 1,
    ],
    [
        144, 137, 130, 124, 113, 107, 101, 95, 88, 82, 76, 70, 64, 57, 51, 45, 39, 33, 26, 15, 1,
    ],
    [
        152, 145, 138, 132, 123, 117, 111, 105, 98, 92, 86, 80, 74, 67, 61, 55, 49, 43, 36, 20, 1,
    ],
    [
        162, 155, 148, 142, 133, 127, 121, 115, 108, 102, 96, 90, 84, 77, 71, 65, 59, 53, 46, 30, 1,
    ],
    [
        172, 165, 158, 152, 143, 137, 131, 125, 118, 112, 106, 100, 94, 87, 81, 75, 69, 63, 56, 45,
        20,
    ],
    [
        200, 200, 200, 200, 200, 200, 200, 200, 198, 193, 188, 183, 178, 173, 168, 163, 158, 153,
        148, 129, 104,
    ],
];

/// RFC 6716 Appendix A celt/quant_bands.c, `eMeans`, extracted by tools/appendix_a_tables.py.
#[rustfmt::skip]
pub const E_MEANS: [f32; 25] = [
    6.4375, 6.25, 5.75, 5.3125, 5.0625, 4.8125,
    4.5, 4.375, 4.875, 4.6875, 4.5625, 4.4375,
    4.875, 4.625, 4.3125, 4.5, 4.375, 4.625,
    4.75, 4.4375, 3.75, 3.75, 3.75, 3.75,
    3.75,
];

/// RFC 6716 Appendix A celt/quant_bands.c, `pred_coef`, extracted by tools/appendix_a_tables.py.
pub const PRED_COEF: [f32; 4] = [0.8984375, 0.796875, 0.6484375, 0.5];

/// RFC 6716 Appendix A celt/quant_bands.c, `beta_coef`, extracted by tools/appendix_a_tables.py.
pub const BETA_COEF: [f32; 4] = [0.9200134, 0.6799927, 0.36999512, 0.2000122];

/// RFC 6716 Appendix A celt/quant_bands.c, `beta_intra`, extracted by tools/appendix_a_tables.py.
pub const BETA_INTRA: f32 = 0.1499939;

/// RFC 6716 Appendix A celt/quant_bands.c, `e_prob_model`, extracted by tools/appendix_a_tables.py.
pub const E_PROB_MODEL: [[[u8; 42]; 2]; 4] = [
    [
        [
            72, 127, 65, 129, 66, 128, 65, 128, 64, 128, 62, 128, 64, 128, 64, 128, 92, 78, 92, 79,
            92, 78, 90, 79, 116, 41, 115, 40, 114, 40, 132, 26, 132, 26, 145, 17, 161, 12, 176, 10,
            177, 11,
        ],
        [
            24, 179, 48, 138, 54, 135, 54, 132, 53, 134, 56, 133, 55, 132, 55, 132, 61, 114, 70,
            96, 74, 88, 75, 88, 87, 74, 89, 66, 91, 67, 100, 59, 108, 50, 120, 40, 122, 37, 97, 43,
            78, 50,
        ],
    ],
    [
        [
            83, 78, 84, 81, 88, 75, 86, 74, 87, 71, 90, 73, 93, 74, 93, 74, 109, 40, 114, 36, 117,
            34, 117, 34, 143, 17, 145, 18, 146, 19, 162, 12, 165, 10, 178, 7, 189, 6, 190, 8, 177,
            9,
        ],
        [
            23, 178, 54, 115, 63, 102, 66, 98, 69, 99, 74, 89, 71, 91, 73, 91, 78, 89, 86, 80, 92,
            66, 93, 64, 102, 59, 103, 60, 104, 60, 117, 52, 123, 44, 138, 35, 133, 31, 97, 38, 77,
            45,
        ],
    ],
    [
        [
            61, 90, 93, 60, 105, 42, 107, 41, 110, 45, 116, 38, 113, 38, 112, 38, 124, 26, 132, 27,
            136, 19, 140, 20, 155, 14, 159, 16, 158, 18, 170, 13, 177, 10, 187, 8, 192, 6, 175, 9,
            159, 10,
        ],
        [
            21, 178, 59, 110, 71, 86, 75, 85, 84, 83, 91, 66, 88, 73, 87, 72, 92, 75, 98, 72, 105,
            58, 107, 54, 115, 52, 114, 55, 112, 56, 129, 51, 132, 40, 150, 33, 140, 29, 98, 35, 77,
            42,
        ],
    ],
    [
        [
            42, 121, 96, 66, 108, 43, 111, 40, 117, 44, 123, 32, 120, 36, 119, 33, 127, 33, 134,
            34, 139, 21, 147, 23, 152, 20, 158, 25, 154, 26, 166, 21, 173, 16, 184, 13, 184, 10,
            150, 13, 139, 15,
        ],
        [
            22, 178, 63, 114, 74, 82, 84, 83, 92, 82, 103, 62, 96, 72, 96, 67, 101, 73, 107, 72,
            113, 55, 118, 52, 125, 52, 118, 52, 117, 55, 135, 49, 137, 39, 157, 32, 145, 29, 97,
            33, 77, 40,
        ],
    ],
];

/// RFC 6716 Appendix A celt/quant_bands.c, `small_energy_icdf`, extracted by tools/appendix_a_tables.py.
pub const SMALL_ENERGY_ICDF: [u8; 3] = [2, 1, 0];

/// RFC 6716 Appendix A celt/celt.c, `trim_icdf`, extracted by tools/appendix_a_tables.py.
pub const TRIM_ICDF: [u8; 11] = [126, 124, 119, 109, 87, 41, 19, 9, 4, 2, 0];

/// RFC 6716 Appendix A celt/celt.c, `spread_icdf`, extracted by tools/appendix_a_tables.py.
pub const SPREAD_ICDF: [u8; 4] = [25, 23, 2, 0];

/// RFC 6716 Appendix A celt/celt.c, `tapset_icdf`, extracted by tools/appendix_a_tables.py.
pub const TAPSET_ICDF: [u8; 3] = [2, 1, 0];

/// RFC 6716 Appendix A celt/celt.c, `tf_select_table`, extracted by tools/appendix_a_tables.py.
pub const TF_SELECT_TABLE: [[i8; 8]; 4] = [
    [0, -1, 0, -1, 0, -1, 0, -1],
    [0, -1, 0, -2, 1, 0, 1, -1],
    [0, -2, 0, -3, 2, 0, 1, -1],
    [0, -2, 0, -3, 3, 0, 1, -1],
];

/// RFC 6716 Appendix A celt/celt.c, `gains`, extracted by tools/appendix_a_tables.py.
pub const COMB_FILTER_GAINS: [[f32; 3]; 3] = [
    [0.30664062, 0.21704102, 0.12963867],
    [0.4638672, 0.2680664, 0.0],
    [0.7998047, 0.100097656, 0.0],
];

/// RFC 6716 Appendix A celt/vq.c, `SPREAD_FACTOR`, extracted by tools/appendix_a_tables.py.
pub const SPREAD_FACTOR: [i32; 3] = [15, 10, 5];

/// RFC 6716 Appendix A celt/static_modes_float.h, `window120`, extracted by tools/appendix_a_tables.py.
#[rustfmt::skip]
pub const WINDOW120: [f32; 120] = [
    6.7286965e-5, 0.00060551346, 0.001681597, 0.0032947962, 0.0054439944, 0.008127692,
    0.011344001, 0.015090633, 0.019364886, 0.024163635, 0.029483315, 0.035319906,
    0.04166891, 0.04852535, 0.055883717, 0.063737996, 0.07208162, 0.08090743,
    0.0902077, 0.09997411, 0.11019769, 0.12086883, 0.13197729, 0.14351214,
    0.15546177, 0.1678139, 0.1805555, 0.1936729, 0.20715171, 0.22097681,
    0.23513243, 0.24960208, 0.2643686, 0.27941418, 0.2947204, 0.3102682,
    0.32603788, 0.3420093, 0.35816178, 0.37447408, 0.39092463, 0.40749142,
    0.42415214, 0.44088423, 0.45766485, 0.47447103, 0.49127978, 0.50806797,
    0.52481264, 0.5414908, 0.5580797, 0.574557, 0.5909005, 0.6070884,
    0.6230995, 0.63891304, 0.65450895, 0.66986775, 0.6849708, 0.6998001,
    0.7143387, 0.7285705, 0.74248046, 0.7560542, 0.76927894, 0.7821426,
    0.7946343, 0.80674446, 0.8184646, 0.8297873, 0.8407067, 0.8512178,
    0.861317, 0.87100184, 0.88027114, 0.8891248, 0.897564, 0.90559095,
    0.913209, 0.9204227, 0.9272374, 0.93365955, 0.93969655, 0.9453567,
    0.9506491, 0.9555835, 0.9601707, 0.9644217, 0.9683485, 0.97196335,
    0.97527903, 0.97830886, 0.98106617, 0.9835648, 0.9858187, 0.9878419,
    0.9896486, 0.9912527, 0.9926685, 0.9939097, 0.99499005, 0.995923,
    0.9967216, 0.99739873, 0.99796665, 0.9984373, 0.998822, 0.99913144,
    0.99937606, 0.99956524, 0.999708, 0.9998125, 0.99988616, 0.9999356,
    0.999967, 0.99998516, 0.9999946, 0.99999857, 0.9999998, 1.0,
];
