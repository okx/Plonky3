//! The Poseidon2 permutation.
//!
//! This implementation was based upon the following resources:
//! - https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2.rs
//! - https://eprint.iacr.org/2023/323.pdf

// #![no_std]

extern crate alloc;

mod diffusion;
mod matrix;
mod round_numbers;
use alloc::vec::Vec;

pub use diffusion::{matmul_internal, DiffusionPermutation};
pub use matrix::*;
use p3_field::{AbstractField, PrimeField, PrimeField64};
use p3_symmetric::{CryptographicPermutation, Permutation};
use rand::distributions::{Distribution, Standard};
use rand::Rng;
pub use round_numbers::poseidon2_round_numbers_128;

const SUPPORTED_WIDTHS: [usize; 8] = [2, 3, 4, 8, 12, 16, 20, 24];

/// The Poseidon2 permutation.
#[derive(Clone, Debug)]
pub struct Poseidon2<F, MdsLight, Diffusion, const WIDTH: usize, const D: u64> {
    /// The number of external rounds.
    rounds_f: usize,

    /// The external round constants.
    external_constants: Vec<[F; WIDTH]>,

    /// The linear layer used in External Rounds. Should be either MDS or a
    /// circulant matrix based off an MDS matrix of size 4.
    external_linear_layer: MdsLight,

    /// The number of internal rounds.
    rounds_p: usize,

    /// The internal round constants.
    internal_constants: Vec<F>,

    /// The linear layer used in internal rounds (only needs diffusion property, not MDS).
    internal_linear_layer: Diffusion,
}

impl<F, MdsLight, Diffusion, const WIDTH: usize, const D: u64>
    Poseidon2<F, MdsLight, Diffusion, WIDTH, D>
where
    F: PrimeField,
{
    /// Create a new Poseidon2 configuration.
    pub fn new(
        rounds_f: usize,
        external_constants: Vec<[F; WIDTH]>,
        external_linear_layer: MdsLight,
        rounds_p: usize,
        internal_constants: Vec<F>,
        internal_linear_layer: Diffusion,
    ) -> Self {
        assert!(SUPPORTED_WIDTHS.contains(&WIDTH));
        Self {
            rounds_f,
            external_constants,
            external_linear_layer,
            rounds_p,
            internal_constants,
            internal_linear_layer,
        }
    }

    /// Create a new Poseidon2 configuration with random parameters.
    pub fn new_from_rng<R: Rng>(
        rounds_f: usize,
        external_linear_layer: MdsLight,
        rounds_p: usize,
        internal_linear_layer: Diffusion,
        rng: &mut R,
    ) -> Self
    where
        Standard: Distribution<F> + Distribution<[F; WIDTH]>,
    {
        let external_constants = rng
            .sample_iter(Standard)
            .take(rounds_f)
            .collect::<Vec<[F; WIDTH]>>();
        let internal_constants = rng.sample_iter(Standard).take(rounds_p).collect::<Vec<F>>();

        Self {
            rounds_f,
            external_constants,
            external_linear_layer,
            rounds_p,
            internal_constants,
            internal_linear_layer,
        }
    }

    #[inline]
    fn add_rc<AF>(&self, state: &mut [AF; WIDTH], rc: &[AF::F; WIDTH])
    where
        AF: AbstractField<F = F>,
    {
        state
            .iter_mut()
            .zip(rc)
            .for_each(|(a, b)| *a += AF::from_f(*b));
    }

    #[inline]
    fn sbox_p<AF>(&self, input: &AF) -> AF
    where
        AF: AbstractField<F = F>,
    {
        input.exp_const_u64::<D>()
    }

    #[inline]
    fn sbox<AF>(&self, state: &mut [AF; WIDTH])
    where
        AF: AbstractField<F = F>,
    {
        state.iter_mut().for_each(|a| *a = self.sbox_p(a));
    }
}

impl<F, MdsLight, Diffusion, const WIDTH: usize, const D: u64>
    Poseidon2<F, MdsLight, Diffusion, WIDTH, D>
where
    F: PrimeField64,
{
    /// Create a new Poseidon2 configuration with 128 bit security and random rounds constants.
    pub fn new_from_rng_128<R: Rng>(
        external_linear_layer: MdsLight,
        internal_linear_layer: Diffusion,
        rng: &mut R,
    ) -> Self
    where
        Standard: Distribution<F> + Distribution<[F; WIDTH]>,
    {
        let (rounds_f, rounds_p) = poseidon2_round_numbers_128::<F>(WIDTH, D);
        println!("rounds_f: {:?}, rounds_p: {:?}", rounds_f, rounds_p);

        let external_constants = rng
            .sample_iter(Standard)
            .take(rounds_f)
            .collect::<Vec<[F; WIDTH]>>();
        let internal_constants = rng.sample_iter(Standard).take(rounds_p).collect::<Vec<F>>();

        Self {
            rounds_f,
            external_constants,
            external_linear_layer,
            rounds_p,
            internal_constants,
            internal_linear_layer,
        }
    }
}

impl<AF, MdsLight, Diffusion, const WIDTH: usize, const D: u64> Permutation<[AF; WIDTH]>
    for Poseidon2<AF::F, MdsLight, Diffusion, WIDTH, D>
where
    AF: AbstractField,
    AF::F: PrimeField,
    MdsLight: MdsLightPermutation<AF, WIDTH>,
    Diffusion: DiffusionPermutation<AF, WIDTH>,
{
    fn permute_mut(&self, state: &mut [AF; WIDTH]) {
        // The initial linear layer.
        self.external_linear_layer.permute_mut(state);

        // The first half of the external rounds.
        let rounds_f_half = self.rounds_f / 2;
        for r in 0..rounds_f_half {
            self.add_rc(state, &self.external_constants[r]);
            self.sbox(state);
            self.external_linear_layer.permute_mut(state);
        }

        // The internal rounds.
        for r in 0..self.rounds_p {
            state[0] += AF::from_f(self.internal_constants[r]);
            state[0] = self.sbox_p(&state[0]);
            self.internal_linear_layer.permute_mut(state);
        }

        // The second half of the external rounds.
        for r in rounds_f_half..self.rounds_f {
            self.add_rc(state, &self.external_constants[r]);
            self.sbox(state);
            self.external_linear_layer.permute_mut(state);
        }
    }
}

impl<AF, MdsLight, Diffusion, const WIDTH: usize, const D: u64>
    CryptographicPermutation<[AF; WIDTH]> for Poseidon2<AF::F, MdsLight, Diffusion, WIDTH, D>
where
    AF: AbstractField,
    AF::F: PrimeField,
    MdsLight: MdsLightPermutation<AF, WIDTH>,
    Diffusion: DiffusionPermutation<AF, WIDTH>,
{
}

/// Field: M31, WIDTH 16, TOTAL_ROUNDS: 30
pub static M31_RC_16_30_U32: [[u32; 16]; 30] = [
    [
        (2110014213),
        (3964964605),
        (2190662774),
        (2732996483),
        (640767983),
        (3403899136),
        (1716033721),
        (1606702601),
        (3759873288),
        (1466015491),
        (1498308946),
        (2844375094),
        (3042463841),
        (1969905919),
        (4109944726),
        (3925048366),
    ],
    [
        (3706859504),
        (759122502),
        (3167665446),
        (1131812921),
        (1080754908),
        (4080114493),
        (893583089),
        (2019677373),
        (3128604556),
        (580640471),
        (3277620260),
        (842931656),
        (548879852),
        (3608554714),
        (3575647916),
        (81826002),
    ],
    [
        (4289086263),
        (1563933798),
        (1440025885),
        (184445025),
        (2598651360),
        (1396647410),
        (1575877922),
        (3303853401),
        (137125468),
        (765010148),
        (633675867),
        (2037803363),
        (2573389828),
        (1895729703),
        (541515871),
        (1783382863),
    ],
    [
        (2641856484),
        (3035743342),
        (3672796326),
        (245668751),
        (2025460432),
        (201609705),
        (286217151),
        (4093475563),
        (2519572182),
        (3080699870),
        (2762001832),
        (1244250808),
        (606038199),
        (3182740831),
        (73007766),
        (2572204153),
    ],
    [
        (1196780786),
        (3447394443),
        (747167305),
        (2968073607),
        (1053214930),
        (1074411832),
        (4016794508),
        (1570312929),
        (113576933),
        (4042581186),
        (3634515733),
        (1032701597),
        (2364839308),
        (3840286918),
        (888378655),
        (2520191583),
    ],
    [
        (36046858),
        (2927525953),
        (3912129105),
        (4004832531),
        (193772436),
        (1590247392),
        (4125818172),
        (2516251696),
        (4050945750),
        (269498914),
        (1973292656),
        (891403491),
        (1845429189),
        (2611996363),
        (2310542653),
        (4071195740),
    ],
    [
        (3505307391),
        (786445290),
        (3815313971),
        (1111591756),
        (4233279834),
        (2775453034),
        (1991257625),
        (2940505809),
        (2751316206),
        (1028870679),
        (1282466273),
        (1059053371),
        (834521354),
        (138721483),
        (3100410803),
        (3843128331),
    ],
    [
        (3878220780),
        (4058162439),
        (1478942487),
        (799012923),
        (496734827),
        (3521261236),
        (755421082),
        (1361409515),
        (392099473),
        (3178453393),
        (4068463721),
        (7935614),
        (4140885645),
        (2150748066),
        (1685210312),
        (3852983224),
    ],
    [
        (2896943075),
        (3087590927),
        (992175959),
        (970216228),
        (3473630090),
        (3899670400),
        (3603388822),
        (2633488197),
        (2479406964),
        (2420952999),
        (1852516800),
        (4253075697),
        (979699862),
        (1163403191),
        (1608599874),
        (3056104448),
    ],
    [
        (3779109343),
        (536205958),
        (4183458361),
        (1649720295),
        (1444912244),
        (3122230878),
        (384301396),
        (4228198516),
        (1662916865),
        (4082161114),
        (2121897314),
        (1706239958),
        (4166959388),
        (1626054781),
        (3005858978),
        (1431907253),
    ],
    [
        (1418914503),
        (1365856753),
        (3942715745),
        (1429155552),
        (3545642795),
        (3772474257),
        (1621094396),
        (2154399145),
        (826697382),
        (1700781391),
        (3539164324),
        (652815039),
        (442484755),
        (2055299391),
        (1064289978),
        (1152335780),
    ],
    [
        (3417648695),
        (186040114),
        (3475580573),
        (2113941250),
        (1779573826),
        (1573808590),
        (3235694804),
        (2922195281),
        (1119462702),
        (3688305521),
        (1849567013),
        (667446787),
        (753897224),
        (1896396780),
        (3143026334),
        (3829603876),
    ],
    [
        (859661334),
        (3898844357),
        (180258337),
        (2321867017),
        (3599002504),
        (2886782421),
        (3038299378),
        (1035366250),
        (2038912197),
        (2920174523),
        (1277696101),
        (2785700290),
        (3806504335),
        (3518858933),
        (654843672),
        (2127120275),
    ],
    [
        (1548195514),
        (2378056027),
        (390914568),
        (1472049779),
        (1552596765),
        (1905886441),
        (1611959354),
        (3653263304),
        (3423946386),
        (340857935),
        (2208879480),
        (139364268),
        (3447281773),
        (3777813707),
        (55640413),
        (4101901741),
    ],
    [
        (104929687),
        (1459980974),
        (1831234737),
        (457139004),
        (2581487628),
        (2112044563),
        (3567013861),
        (2792004347),
        (576325418),
        (41126132),
        (2713562324),
        (151213722),
        (2891185935),
        (546846420),
        (2939794919),
        (2543469905),
    ],
    [
        (2191909784),
        (3315138460),
        (530414574),
        (1242280418),
        (1211740715),
        (3993672165),
        (2505083323),
        (3845798801),
        (538768466),
        (2063567560),
        (3366148274),
        (1449831887),
        (2408012466),
        (294726285),
        (3943435493),
        (924016661),
    ],
    [
        (3633138367),
        (3222789372),
        (809116305),
        (30100013),
        (2655172876),
        (2564247117),
        (2478649732),
        (4113689151),
        (4120146082),
        (2512308515),
        (650406041),
        (4240012393),
        (2683508708),
        (951073977),
        (3460081988),
        (339124269),
    ],
    [
        (130182653),
        (2755946749),
        (542600513),
        (2816103022),
        (1931786340),
        (2044470840),
        (1709908013),
        (2938369043),
        (3640399693),
        (1374470239),
        (2191149676),
        (2637495682),
        (4236394040),
        (2289358846),
        (3833368530),
        (974546524),
    ],
    [
        (3306659113),
        (2234814261),
        (1188782305),
        (223782844),
        (2248980567),
        (2309786141),
        (2023401627),
        (3278877413),
        (2022138149),
        (575851471),
        (1612560780),
        (3926656936),
        (3318548977),
        (2591863678),
        (188109355),
        (4217723909),
    ],
    [
        (1564209905),
        (2154197895),
        (2459687029),
        (2870634489),
        (1375012945),
        (1529454825),
        (306140690),
        (2855578299),
        (1246997295),
        (3024298763),
        (1915270363),
        (1218245412),
        (2479314020),
        (2989827755),
        (814378556),
        (4039775921),
    ],
    [
        (1165280628),
        (1203983801),
        (3814740033),
        (1919627044),
        (600240215),
        (773269071),
        (486685186),
        (4254048810),
        (1415023565),
        (502840102),
        (4225648358),
        (510217063),
        (166444818),
        (1430745893),
        (1376516190),
        (1775891321),
    ],
    [
        (1170945922),
        (1105391877),
        (261536467),
        (1401687994),
        (1022529847),
        (2476446456),
        (2603844878),
        (3706336043),
        (3463053714),
        (1509644517),
        (588552318),
        (65252581),
        (3696502656),
        (2183330763),
        (3664021233),
        (1643809916),
    ],
    [
        (2922875898),
        (3740690643),
        (3932461140),
        (161156271),
        (2619943483),
        (4077039509),
        (2921201703),
        (2085619718),
        (2065264646),
        (2615693812),
        (3116555433),
        (246100007),
        (4281387154),
        (4046141001),
        (4027749321),
        (111611860),
    ],
    [
        (2066954820),
        (2502099969),
        (2915053115),
        (2362518586),
        (366091708),
        (2083204932),
        (4138385632),
        (3195157567),
        (1318086382),
        (521723799),
        (702443405),
        (2507670985),
        (1760347557),
        (2631999893),
        (1672737554),
        (1060867760),
    ],
    [
        (2359801781),
        (2800231467),
        (3010357035),
        (1035997899),
        (1210110952),
        (1018506770),
        (2799468177),
        (1479380761),
        (1536021911),
        (358993854),
        (579904113),
        (3432144800),
        (3625515809),
        (199241497),
        (4058304109),
        (2590164234),
    ],
    [
        (1688530738),
        (1580733335),
        (2443981517),
        (2206270565),
        (2780074229),
        (2628739677),
        (2940123659),
        (4145206827),
        (3572278009),
        (2779607509),
        (1098718697),
        (1424913749),
        (2224415875),
        (1108922178),
        (3646272562),
        (3935186184),
    ],
    [
        (820046587),
        (1393386250),
        (2665818575),
        (2231782019),
        (672377010),
        (1920315467),
        (1913164407),
        (2029526876),
        (2629271820),
        (384320012),
        (4112320585),
        (3131824773),
        (2347818197),
        (2220997386),
        (1772368609),
        (2579960095),
    ],
    [
        (3544930873),
        (225847443),
        (3070082278),
        (95643305),
        (3438572042),
        (3312856509),
        (615850007),
        (1863868773),
        (803582265),
        (3461976859),
        (2903025799),
        (1482092434),
        (3902972499),
        (3872341868),
        (1530411808),
        (2214923584),
    ],
    [
        (3118792481),
        (2241076515),
        (3983669831),
        (3180915147),
        (3838626501),
        (1921630011),
        (3415351771),
        (2249953859),
        (3755081630),
        (486327260),
        (1227575720),
        (3643869379),
        (2982026073),
        (2466043731),
        (1982634375),
        (3769609014),
    ],
    [
        (2195455495),
        (2596863283),
        (4244994973),
        (1983609348),
        (4019674395),
        (3469982031),
        (1458697570),
        (1593516217),
        (1963896497),
        (3115309118),
        (1659132465),
        (2536770756),
        (3059294171),
        (2618031334),
        (2040903247),
        (3799795076),
    ],
];
