/*
	conversions
*/
#ifndef HAVE_EXPLICIT_BZERO
# define explicit_bzero(b, len) memset((b), 0, (len))
#endif

DONNA_INLINE static void
ge25519_p1p1_to_partial(ge25519 *r, const ge25519_p1p1 *p) {
	curve25519_mul(r->x, p->x, p->t);
	curve25519_mul(r->y, p->y, p->z);
	curve25519_mul(r->z, p->z, p->t); 
}

DONNA_INLINE static void
ge25519_p1p1_to_full(ge25519 *r, const ge25519_p1p1 *p) {
	curve25519_mul(r->x, p->x, p->t);
	curve25519_mul(r->y, p->y, p->z);
	curve25519_mul(r->z, p->z, p->t); 
	curve25519_mul(r->t, p->x, p->y); 
}

static void
ge25519_full_to_pniels(ge25519_pniels *p, const ge25519 *r) {
	curve25519_sub(p->ysubx, r->y, r->x);
	curve25519_add(p->xaddy, r->y, r->x);
	curve25519_copy(p->z, r->z);
	curve25519_mul(p->t2d, r->t, ge25519_ec2d);
}

/*
	adding & doubling
*/

static void
ge25519_add_p1p1(ge25519_p1p1 *r, const ge25519 *p, const ge25519 *q) {
	bignum25519 a,b,c,d,t,u;

	curve25519_sub(a, p->y, p->x);
	curve25519_add(b, p->y, p->x);
	curve25519_sub(t, q->y, q->x);
	curve25519_add(u, q->y, q->x);
	curve25519_mul(a, a, t);
	curve25519_mul(b, b, u);
	curve25519_mul(c, p->t, q->t);
	curve25519_mul(c, c, ge25519_ec2d);
	curve25519_mul(d, p->z, q->z);
	curve25519_add(d, d, d);
	curve25519_sub(r->x, b, a);
	curve25519_add(r->y, b, a);
	curve25519_add_after_basic(r->z, d, c);
	curve25519_sub_after_basic(r->t, d, c);
}


static void
ge25519_double_p1p1(ge25519_p1p1 *r, const ge25519 *p) {
	bignum25519 a,b,c;

	curve25519_square(a, p->x);
	curve25519_square(b, p->y);
	curve25519_square(c, p->z);
	curve25519_add_reduce(c, c, c);
	curve25519_add(r->x, p->x, p->y);
	curve25519_square(r->x, r->x);
	curve25519_add(r->y, b, a);
	curve25519_sub(r->z, b, a);
	curve25519_sub_after_basic(r->x, r->x, r->y);
	curve25519_sub_after_basic(r->t, c, r->z);
}

static void
ge25519_nielsadd2_p1p1(ge25519_p1p1 *r, const ge25519 *p, const ge25519_niels *q, unsigned char signbit) {
	const bignum25519 *qb = (const bignum25519 *)q;
	bignum25519 *rb = (bignum25519 *)r;
	bignum25519 a,b,c;

	curve25519_sub(a, p->y, p->x);
	curve25519_add(b, p->y, p->x);
	curve25519_mul(a, a, qb[signbit]); /* x for +, y for - */
	curve25519_mul(r->x, b, qb[signbit^1]); /* y for +, x for - */
	curve25519_add(r->y, r->x, a);
	curve25519_sub(r->x, r->x, a);
	curve25519_mul(c, p->t, q->t2d);
	curve25519_add_reduce(r->t, p->z, p->z);
	curve25519_copy(r->z, r->t);
	curve25519_add(rb[2+signbit], rb[2+signbit], c); /* z for +, t for - */
	curve25519_sub(rb[2+(signbit^1)], rb[2+(signbit^1)], c); /* t for +, z for - */
}

static void
ge25519_pnielsadd_p1p1(ge25519_p1p1 *r, const ge25519 *p, const ge25519_pniels *q, unsigned char signbit) {
	const bignum25519 *qb = (const bignum25519 *)q;
	bignum25519 *rb = (bignum25519 *)r;
	bignum25519 a,b,c;

	curve25519_sub(a, p->y, p->x);
	curve25519_add(b, p->y, p->x);
	curve25519_mul(a, a, qb[signbit]); /* ysubx for +, xaddy for - */
	curve25519_mul(r->x, b, qb[signbit^1]); /* xaddy for +, ysubx for - */
	curve25519_add(r->y, r->x, a);
	curve25519_sub(r->x, r->x, a);
	curve25519_mul(c, p->t, q->t2d);
	curve25519_mul(r->t, p->z, q->z);
	curve25519_add_reduce(r->t, r->t, r->t);
	curve25519_copy(r->z, r->t);
	curve25519_add(rb[2+signbit], rb[2+signbit], c); /* z for +, t for - */
	curve25519_sub(rb[2+(signbit^1)], rb[2+(signbit^1)], c); /* t for +, z for - */
}

static void
ge25519_double_partial(ge25519 *r, const ge25519 *p) {
	ge25519_p1p1 t;
	ge25519_double_p1p1(&t, p);
	ge25519_p1p1_to_partial(r, &t);
}

static void
ge25519_double(ge25519 *r, const ge25519 *p) {
	ge25519_p1p1 t;
	ge25519_double_p1p1(&t, p);
	ge25519_p1p1_to_full(r, &t);
}

STATIC void
ge25519_add(ge25519 *r, const ge25519 *p,  const ge25519 *q) {
	ge25519_p1p1 t;
	ge25519_add_p1p1(&t, p, q);
	ge25519_p1p1_to_full(r, &t);
}

static void
ge25519_nielsadd2(ge25519 *r, const ge25519_niels *q) {
	bignum25519 a,b,c,e,f,g,h;

	curve25519_sub(a, r->y, r->x);
	curve25519_add(b, r->y, r->x);
	curve25519_mul(a, a, q->ysubx);
	curve25519_mul(e, b, q->xaddy);
	curve25519_add(h, e, a);
	curve25519_sub(e, e, a);
	curve25519_mul(c, r->t, q->t2d);
	curve25519_add(f, r->z, r->z);
	curve25519_add_after_basic(g, f, c);
	curve25519_sub_after_basic(f, f, c);
	curve25519_mul(r->x, e, f);
	curve25519_mul(r->y, h, g);
	curve25519_mul(r->z, g, f);
	curve25519_mul(r->t, e, h);
}

static void
ge25519_pnielsadd(ge25519_pniels *r, const ge25519 *p, const ge25519_pniels *q) {
	bignum25519 a,b,c,x,y,z,t;

	curve25519_sub(a, p->y, p->x);
	curve25519_add(b, p->y, p->x);
	curve25519_mul(a, a, q->ysubx);
	curve25519_mul(x, b, q->xaddy);
	curve25519_add(y, x, a);
	curve25519_sub(x, x, a);
	curve25519_mul(c, p->t, q->t2d);
	curve25519_mul(t, p->z, q->z);
	curve25519_add(t, t, t);
	curve25519_add_after_basic(z, t, c);
	curve25519_sub_after_basic(t, t, c);
	curve25519_mul(r->xaddy, x, t);
	curve25519_mul(r->ysubx, y, z);
	curve25519_mul(r->z, z, t);
	curve25519_mul(r->t2d, x, y);
	curve25519_copy(y, r->ysubx);
	curve25519_sub(r->ysubx, r->ysubx, r->xaddy);
	curve25519_add(r->xaddy, r->xaddy, y);
	curve25519_mul(r->t2d, r->t2d, ge25519_ec2d);
}


/*
	pack & unpack
*/

STATIC void
ge25519_pack(unsigned char r[32], const ge25519 *p) {
	bignum25519 tx, ty, zi;
	unsigned char parity[32];
	curve25519_recip(zi, p->z);
	curve25519_mul(tx, p->x, zi);
	curve25519_mul(ty, p->y, zi);
	curve25519_contract(r, ty);
	curve25519_contract(parity, tx);
	r[31] ^= ((parity[0] & 1) << 7);
}

STATIC int ge25519_unpack_negative_vartime(ge25519 *r, const unsigned char p[32]) {
	static const unsigned long zero[32/sizeof(long)] = {0L};
	static const bignum25519 one = {1};
	unsigned char parity = p[31] >> 7;
	unsigned long check[32/sizeof(long)];
	bignum25519 t, root, num, den, d3;

	curve25519_expand(r->y, p);
	curve25519_copy(r->z, one);
	curve25519_square(num, r->y); /* x = y^2 */
	curve25519_mul(den, num, ge25519_ecd); /* den = dy^2 */
	curve25519_sub_reduce(num, num, r->z); /* x = y^1 - 1 */
	curve25519_add(den, den, r->z); /* den = dy^2 + 1 */

	/* Computation of sqrt(num/den) */
	/* 1.: computation of num^((p-5)/8)*den^((7p-35)/8) = (num*den^7)^((p-5)/8) */
	curve25519_square(t, den);
	curve25519_mul(d3, t, den);
	curve25519_square(r->x, d3);
	curve25519_mul(r->x, r->x, den);
	curve25519_mul(r->x, r->x, num);
	curve25519_pow_two252m3(r->x, r->x);

	/* 2. computation of r->x = num * den^3 * (num*den^7)^((p-5)/8) */
	curve25519_mul(r->x, r->x, d3);
	curve25519_mul(r->x, r->x, num);

	/* 3. Check if either of the roots works: */
	curve25519_square(t, r->x);
	curve25519_mul(t, t, den);
	curve25519_sub_reduce(root, t, num);
	curve25519_contract((unsigned char*)check, root);
	if (!ed25519_verify(check, zero, 32)) {
		curve25519_add_reduce(t, t, num);
		curve25519_contract((unsigned char*)check, t);
		if (!ed25519_verify(check, zero, 32))
			return 0;
		curve25519_mul(r->x, r->x, ge25519_sqrtneg1);
	}

	curve25519_contract((unsigned char*)check, r->x);
	if ((check[0] & 1) == parity) {
		curve25519_copy(t, r->x);
		curve25519_neg(r->x, t);
	}
	curve25519_mul(r->t, r->x, r->y);
	return 1;
}


/*
	scalarmults
*/

DONNA_INLINE static void ge25519_set_neutral(ge25519 *r)
{
 	memset(r, 0, sizeof(ge25519));
	r->y[0] = 1;
	r->z[0] = 1;
}

#define S1_SWINDOWSIZE 5
#define S1_TABLE_SIZE (1<<(S1_SWINDOWSIZE-2))
#define S2_SWINDOWSIZE 7
#define S2_TABLE_SIZE (1<<(S2_SWINDOWSIZE-2))

/* computes [s1]p1 + [s2]base */
STATIC void ge25519_double_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const bignum256modm s1, const bignum256modm s2) {
	signed char slide1[256], slide2[256];
	ge25519_pniels pre1[S1_TABLE_SIZE];
#if defined(ED25519_64BIT)
	const ge25519_niels table2[32] = {
	  {{0x00003905d740913e,0x0000ba2817d673a2,0x00023e2827f4e67c,0x000133d2e0c21a34,0x00044fd2f9298f81},{0x000493c6f58c3b85,0x0000df7181c325f7,0x0000f50b0b3e4cb7,0x0005329385a44c32,0x00007cf9d3a33d4b},{0x00011205877aaa68,0x000479955893d579,0x00050d66309b67a0,0x0002d42d0dbee5ee,0x0006f117b689f0c6}},
	  {{0x00011fe8a4fcd265,0x0007bcb8374faacc,0x00052f5af4ef4d4f,0x0005314098f98d10,0x0002ab91587555bd},{0x0005b0a84cee9730,0x00061d10c97155e4,0x0004059cc8096a10,0x00047a608da8014f,0x0007a164e1b9a80f},{0x0006933f0dd0d889,0x00044386bb4c4295,0x0003cb6d3162508c,0x00026368b872a2c6,0x0005a2826af12b9b}},
	  {{0x000182c3a447d6ba,0x00022964e536eff2,0x000192821f540053,0x0002f9f19e788e5c,0x000154a7e73eb1b5},{0x0002bc4408a5bb33,0x000078ebdda05442,0x0002ffb112354123,0x000375ee8df5862d,0x0002945ccf146e20},{0x0003dbf1812a8285,0x0000fa17ba3f9797,0x0006f69cb49c3820,0x00034d5a0db3858d,0x00043aabe696b3bb}},
	  {{0x00072c9aaa3221b1,0x000267774474f74d,0x000064b0e9b28085,0x0003f04ef53b27c9,0x0001d6edd5d2e531},{0x00025cd0944ea3bf,0x00075673b81a4d63,0x000150b925d1c0d4,0x00013f38d9294114,0x000461bea69283c9},{0x00036dc801b8b3a2,0x0000e0a7d4935e30,0x0001deb7cecc0d7d,0x000053a94e20dd2c,0x0007a9fbb1c6a0f9}},
	  {{0x0006217e039d8064,0x0006dea408337e6d,0x00057ac112628206,0x000647cb65e30473,0x00049c05a51fadc9},{0x0006678aa6a8632f,0x0005ea3788d8b365,0x00021bd6d6994279,0x0007ace75919e4e3,0x00034b9ed338add7},{0x0004e8bf9045af1b,0x000514e33a45e0d6,0x0007533c5b8bfe0f,0x000583557b7e14c9,0x00073c172021b008}},
	  {{0x00075b0249864348,0x00052ee11070262b,0x000237ae54fb5acd,0x0003bfd1d03aaab5,0x00018ab598029d5c},{0x000700848a802ade,0x0001e04605c4e5f7,0x0005c0d01b9767fb,0x0007d7889f42388b,0x0004275aae2546d8},{0x00032cc5fd6089e9,0x000426505c949b05,0x00046a18880c7ad2,0x0004a4221888ccda,0x0003dc65522b53df}},
	  {{0x0007013b327fbf93,0x0001336eeded6a0d,0x0002b565a2bbf3af,0x000253ce89591955,0x0000267882d17602},{0x0000c222a2007f6d,0x000356b79bdb77ee,0x00041ee81efe12ce,0x000120a9bd07097d,0x000234fd7eec346f},{0x0000a119732ea378,0x00063bf1ba8e2a6c,0x00069f94cc90df9a,0x000431d1779bfc48,0x000497ba6fdaa097}},
	  {{0x0003cd86468ccf0b,0x00048553221ac081,0x0006c9464b4e0a6e,0x00075fba84180403,0x00043b5cd4218d05},{0x0006cc0313cfeaa0,0x0001a313848da499,0x0007cb534219230a,0x00039596dedefd60,0x00061e22917f12de},{0x0002762f9bd0b516,0x0001c6e7fbddcbb3,0x00075909c3ace2bd,0x00042101972d3ec9,0x000511d61210ae4d}},
	  {{0x000386484420de87,0x0002d6b25db68102,0x000650b4962873c0,0x0004081cfd271394,0x00071a7fe6fe2482},{0x000676ef950e9d81,0x0001b81ae089f258,0x00063c4922951883,0x0002f1d54d9b3237,0x0006d325924ddb85},{0x000182b8a5c8c854,0x00073fcbe5406d8e,0x0005de3430cff451,0x000554b967ac8c41,0x0004746c4b6559ee}},
	  {{0x000546c864741147,0x0003a1df99092690,0x0001ca8cc9f4d6bb,0x00036b7fc9cd3b03,0x000219663497db5e},{0x00077b3c6dc69a2b,0x0004edf13ec2fa6e,0x0004e85ad77beac8,0x0007dba2b28e7bda,0x0005c9a51de34fe9},{0x0000f1cf79f10e67,0x00043ccb0a2b7ea2,0x00005089dfff776a,0x0001dd84e1d38b88,0x0004804503c60822}},
	  {{0x000021d23a36d175,0x0004fd3373c6476d,0x00020e291eeed02a,0x00062f2ecf2e7210,0x000771e098858de4},{0x00049ed02ca37fc7,0x000474c2b5957884,0x0005b8388e816683,0x0004b6c454b76be4,0x000553398a516506},{0x0002f5d278451edf,0x000730b133997342,0x0006965420eb6975,0x000308a3bfa516cf,0x0005a5ed1d68ff5a}},
	  {{0x0005e0c558527359,0x0003395b73afd75c,0x000072afa4e4b970,0x00062214329e0f6d,0x000019b60135fefd},{0x0005122afe150e83,0x0004afc966bb0232,0x0001c478833c8268,0x00017839c3fc148f,0x00044acb897d8bf9},{0x000068145e134b83,0x0001e4860982c3cc,0x000068fb5f13d799,0x0007c9283744547e,0x000150c49fde6ad2}},
	  {{0x0001863c9cdca868,0x0003770e295a1709,0x0000d85a3720fd13,0x0005e0ff1f71ab06,0x00078a6d7791e05f},{0x0003f29509471138,0x000729eeb4ca31cf,0x00069c22b575bfbc,0x0004910857bce212,0x0006b2b5a075bb99},{0x0007704b47a0b976,0x0002ae82e91aab17,0x00050bd6429806cd,0x00068055158fd8ea,0x000725c7ffc4ad55}},
	  {{0x00002bf71cd098c0,0x00049dabcc6cd230,0x00040a6533f905b2,0x000573efac2eb8a4,0x0004cd54625f855f},{0x00026715d1cf99b2,0x0002205441a69c88,0x000448427dcd4b54,0x0001d191e88abdc5,0x000794cc9277cb1f},{0x0006c426c2ac5053,0x0005a65ece4b095e,0x0000c44086f26bb6,0x0007429568197885,0x0007008357b6fcc8}},
	  {{0x00039fbb82584a34,0x00047a568f257a03,0x00014d88091ead91,0x0002145b18b1ce24,0x00013a92a3669d6d},{0x0000672738773f01,0x000752bf799f6171,0x0006b4a6dae33323,0x0007b54696ead1dc,0x00006ef7e9851ad0},{0x0003771cc0577de5,0x0003ca06bb8b9952,0x00000b81c5d50390,0x00043512340780ec,0x0003c296ddf8a2af}},
	  {{0x00034d2ebb1f2541,0x0000e815b723ff9d,0x000286b416e25443,0x0000bdfe38d1bee8,0x0000a892c7007477},{0x000515f9d914a713,0x00073191ff2255d5,0x00054f5cc2a4bdef,0x0003dd57fc118bcf,0x0007a99d393490c7},{0x0002ed2436bda3e8,0x00002afd00f291ea,0x0000be7381dea321,0x0003e952d4b2b193,0x000286762d28302f}},
	  {{0x00058e2bce2ef5bd,0x00068ce8f78c6f8a,0x0006ee26e39261b2,0x00033d0aa50bcf9d,0x0007686f2a3d6f17},{0x000036093ce35b25,0x0003b64d7552e9cf,0x00071ee0fe0b8460,0x00069d0660c969e5,0x00032f1da046a9d9},{0x000512a66d597c6a,0x0000609a70a57551,0x000026c08a3c464c,0x0004531fc8ee39e1,0x000561305f8a9ad2}},
	  {{0x0002cc28e7b0c0d5,0x00077b60eb8a6ce4,0x0004042985c277a6,0x000636657b46d3eb,0x000030a1aef2c57c},{0x0004978dec92aed1,0x000069adae7ca201,0x00011ee923290f55,0x00069641898d916c,0x00000aaec53e35d4},{0x0001f773003ad2aa,0x000005642cc10f76,0x00003b48f82cfca6,0x0002403c10ee4329,0x00020be9c1c24065}},
	  {{0x0000e44ae2025e60,0x0005f97b9727041c,0x0005683472c0ecec,0x000188882eb1ce7c,0x00069764c545067e},{0x000387d8249673a6,0x0005bea8dc927c2a,0x0005bd8ed5650ef0,0x0000ef0e3fcd40e1,0x000750ab3361f0ac},{0x00023283a2f81037,0x000477aff97e23d1,0x0000b8958dbcbb68,0x0000205b97e8add6,0x00054f96b3fb7075}},
	  {{0x0005afc616b11ecd,0x00039f4aec8f22ef,0x0003b39e1625d92e,0x0005f85bd4508873,0x00078e6839fbe85d},{0x0005f20429669279,0x00008fafae4941f5,0x00015d83c4eb7688,0x0001cf379eca4146,0x0003d7fe9c52bb75},{0x00032df737b8856b,0x0000608342f14e06,0x0003967889d74175,0x0001211907fba550,0x00070f268f350088}},
	  {{0x0004112070dcf355,0x0007dcff9c22e464,0x00054ada60e03325,0x00025cd98eef769a,0x000404e56c039b8c},{0x00064583b1805f47,0x00022c1baf832cd0,0x000132c01bd4d717,0x0004ecf4c3a75b8f,0x0007c0d345cfad88},{0x00071f4b8c78338a,0x00062cfc16bc2b23,0x00017cf51280d9aa,0x0003bbae5e20a95a,0x00020d754762aaec}},
	  {{0x0004feb135b9f543,0x00063bd192ad93ae,0x00044e2ea612cdf7,0x000670f4991583ab,0x00038b8ada8790b4},{0x0007c36fc73bb758,0x0004a6c797734bd1,0x0000ef248ab3950e,0x00063154c9a53ec8,0x0002b8f1e46f3cee},{0x00004a9cdf51f95d,0x0005d963fbd596b8,0x00022d9b68ace54a,0x0004a98e8836c599,0x000049aeb32ceba1}},
	  {{0x00067d3c63dcfe7e,0x000112f0adc81aee,0x00053df04c827165,0x0002fe5b33b430f0,0x00051c665e0c8d62},{0x00007d0b75fc7931,0x00016f4ce4ba754a,0x0005ace4c03fbe49,0x00027e0ec12a159c,0x000795ee17530f67},{0x00025b0a52ecbd81,0x0005dc0695fce4a9,0x0003b928c575047d,0x00023bf3512686e5,0x0006cd19bf49dc54}},
	  {{0x0007619052179ca3,0x0000c16593f0afd0,0x000265c4795c7428,0x00031c40515d5442,0x0007520f3db40b2e},{0x0006612165afc386,0x0001171aa36203ff,0x0002642ea820a8aa,0x0001f3bb7b313f10,0x0005e01b3a7429e4},{0x00050be3d39357a1,0x0003ab33d294a7b6,0x0004c479ba59edb3,0x0004c30d184d326f,0x00071092c9ccef3c}},
	  {{0x0000523f0364918c,0x000687f56d638a7b,0x00020796928ad013,0x0005d38405a54f33,0x0000ea15b03d0257},{0x0003d8ac74051dcf,0x00010ab6f543d0ad,0x0005d0f3ac0fda90,0x0005ef1d2573e5e4,0x0004173a5bb7137a},{0x00056e31f0f9218a,0x0005635f88e102f8,0x0002cbc5d969a5b8,0x000533fbc98b347a,0x0005fc565614a4e3}},
	  {{0x0006570dc46d7ae5,0x00018a9f1b91e26d,0x000436b6183f42ab,0x000550acaa4f8198,0x00062711c414c454},{0x0002e1e67790988e,0x0001e38b9ae44912,0x000648fbb4075654,0x00028df1d840cd72,0x0003214c7409d466},{0x0001827406651770,0x0004d144f286c265,0x00017488f0ee9281,0x00019e6cdb5c760c,0x0005bea94073ecb8}},
	  {{0x0005bf0912c89be4,0x00062fadcaf38c83,0x00025ec196b3ce2c,0x00077655ff4f017b,0x0003aacd5c148f61},{0x0000ce63f343d2f8,0x0001e0a87d1e368e,0x000045edbc019eea,0x0006979aed28d0d1,0x0004ad0785944f1b},{0x00063b34c3318301,0x0000e0e62d04d0b1,0x000676a233726701,0x00029e9a042d9769,0x0003aff0cb1d9028}},
	  {{0x0005c7eb3a20405e,0x0005fdb5aad930f8,0x0004a757e63b8c47,0x00028e9492972456,0x000110e7e86f4cd2},{0x0006430bf4c53505,0x000264c3e4507244,0x00074c9f19a39270,0x00073f84f799bc47,0x0002ccf9f732bd99},{0x0000d89ed603f5e4,0x00051e1604018af8,0x0000b8eedc4a2218,0x00051ba98b9384d0,0x00005c557e0b9693}},
	  {{0x0001ce311fc97e6f,0x0006023f3fb5db1f,0x0007b49775e8fc98,0x0003ad70adbf5045,0x0006e154c178fe98},{0x0006bbb089c20eb0,0x0006df41fb0b9eee,0x00051087ed87e16f,0x000102db5c9fa731,0x000289fef0841861},{0x00016336fed69abf,0x0004f066b929f9ec,0x0004e9ff9e6c5b93,0x00018c89bc4bb2ba,0x0006afbf642a95ca}},
	  {{0x0000de0c62f5d2c1,0x00049601cf734fb5,0x0006b5c38263f0f6,0x0004623ef5b56d06,0x0000db4b851b9503},{0x00055070f913a8cc,0x000765619eac2bbc,0x0003ab5225f47459,0x00076ced14ab5b48,0x00012c093cedb801},{0x00047f9308b8190f,0x000414235c621f82,0x00031f5ff41a5a76,0x0006736773aab96d,0x00033aa8799c6635}},
	  {{0x0007f51ebd085cf2,0x00012cfa67e3f5e1,0x0001800cf1e3d46a,0x00054337615ff0a8,0x000233c6f29e8e21},{0x0000f588fc156cb1,0x000363414da4f069,0x0007296ad9b68aea,0x0004d3711316ae43,0x000212cd0c1c8d58},{0x0004d5107f18c781,0x00064a4fd3a51a5e,0x0004f4cd0448bb37,0x000671d38543151e,0x0001db7778911914}},
	  {{0x000352397c6bc26f,0x00018a7aa0227bbe,0x0005e68cc1ea5f8b,0x0006fe3e3a7a1d5f,0x00031ad97ad26e2a},{0x00014769dd701ab6,0x00028339f1b4b667,0x0004ab214b8ae37b,0x00025f0aefa0b0fe,0x0007ae2ca8a017d2},{0x000017ed0920b962,0x000187e33b53b6fd,0x00055829907a1463,0x000641f248e0a792,0x0001ed1fc53a6622}}
	};
#else
	const ge25519_niels ALIGN(16) table2[32] = {
	  {{0x0340913e,0x000e4175,0x03d673a2,0x002e8a05,0x03f4e67c,0x008f8a09,0x00c21a34,0x004cf4b8,0x01298f81,0x0113f4be},{0x018c3b85,0x0124f1bd,0x01c325f7,0x0037dc60,0x033e4cb7,0x003d42c2,0x01a44c32,0x014ca4e1,0x03a33d4b,0x001f3e74},{0x037aaa68,0x00448161,0x0093d579,0x011e6556,0x009b67a0,0x0143598c,0x01bee5ee,0x00b50b43,0x0289f0c6,0x01bc45ed}},
	  {{0x00fcd265,0x0047fa29,0x034faacc,0x01ef2e0d,0x00ef4d4f,0x014bd6bd,0x00f98d10,0x014c5026,0x007555bd,0x00aae456},{0x00ee9730,0x016c2a13,0x017155e4,0x01874432,0x00096a10,0x01016732,0x01a8014f,0x011e9823,0x01b9a80f,0x01e85938},{0x01d0d889,0x01a4cfc3,0x034c4295,0x0110e1ae,0x0162508c,0x00f2db4c,0x0072a2c6,0x0098da2e,0x02f12b9b,0x0168a09a}},
	  {{0x0047d6ba,0x0060b0e9,0x0136eff2,0x008a5939,0x03540053,0x0064a087,0x02788e5c,0x00be7c67,0x033eb1b5,0x005529f9},{0x00a5bb33,0x00af1102,0x01a05442,0x001e3af7,0x02354123,0x00bfec44,0x01f5862d,0x00dd7ba3,0x03146e20,0x00a51733},{0x012a8285,0x00f6fc60,0x023f9797,0x003e85ee,0x009c3820,0x01bda72d,0x01b3858d,0x00d35683,0x0296b3bb,0x010eaaf9}},
	  {{0x023221b1,0x01cb26aa,0x0074f74d,0x0099ddd1,0x01b28085,0x00192c3a,0x013b27c9,0x00fc13bd,0x01d2e531,0x0075bb75},{0x004ea3bf,0x00973425,0x001a4d63,0x01d59cee,0x01d1c0d4,0x00542e49,0x01294114,0x004fce36,0x029283c9,0x01186fa9},{0x01b8b3a2,0x00db7200,0x00935e30,0x003829f5,0x02cc0d7d,0x0077adf3,0x0220dd2c,0x0014ea53,0x01c6a0f9,0x01ea7eec}},
	  {{0x039d8064,0x01885f80,0x00337e6d,0x01b7a902,0x02628206,0x015eb044,0x01e30473,0x0191f2d9,0x011fadc9,0x01270169},{0x02a8632f,0x0199e2a9,0x00d8b365,0x017a8de2,0x02994279,0x0086f5b5,0x0119e4e3,0x01eb39d6,0x0338add7,0x00d2e7b4},{0x0045af1b,0x013a2fe4,0x0245e0d6,0x014538ce,0x038bfe0f,0x01d4cf16,0x037e14c9,0x0160d55e,0x0021b008,0x01cf05c8}},
	  {{0x01864348,0x01d6c092,0x0070262b,0x014bb844,0x00fb5acd,0x008deb95,0x003aaab5,0x00eff474,0x00029d5c,0x0062ad66},{0x02802ade,0x01c02122,0x01c4e5f7,0x00781181,0x039767fb,0x01703406,0x0342388b,0x01f5e227,0x022546d8,0x0109d6ab},{0x016089e9,0x00cb317f,0x00949b05,0x01099417,0x000c7ad2,0x011a8622,0x0088ccda,0x01290886,0x022b53df,0x00f71954}},
	  {{0x027fbf93,0x01c04ecc,0x01ed6a0d,0x004cdbbb,0x02bbf3af,0x00ad5968,0x01591955,0x0094f3a2,0x02d17602,0x00099e20},{0x02007f6d,0x003088a8,0x03db77ee,0x00d5ade6,0x02fe12ce,0x0107ba07,0x0107097d,0x00482a6f,0x02ec346f,0x008d3f5f},{0x032ea378,0x0028465c,0x028e2a6c,0x018efc6e,0x0090df9a,0x01a7e533,0x039bfc48,0x010c745d,0x03daa097,0x0125ee9b}},
	  {{0x028ccf0b,0x00f36191,0x021ac081,0x012154c8,0x034e0a6e,0x01b25192,0x00180403,0x01d7eea1,0x00218d05,0x010ed735},{0x03cfeaa0,0x01b300c4,0x008da499,0x0068c4e1,0x0219230a,0x01f2d4d0,0x02defd60,0x00e565b7,0x017f12de,0x018788a4},{0x03d0b516,0x009d8be6,0x03ddcbb3,0x0071b9fe,0x03ace2bd,0x01d64270,0x032d3ec9,0x01084065,0x0210ae4d,0x01447584}},
	  {{0x0020de87,0x00e19211,0x01b68102,0x00b5ac97,0x022873c0,0x01942d25,0x01271394,0x0102073f,0x02fe2482,0x01c69ff9},{0x010e9d81,0x019dbbe5,0x0089f258,0x006e06b8,0x02951883,0x018f1248,0x019b3237,0x00bc7553,0x024ddb85,0x01b4c964},{0x01c8c854,0x0060ae29,0x01406d8e,0x01cff2f9,0x00cff451,0x01778d0c,0x03ac8c41,0x01552e59,0x036559ee,0x011d1b12}},
	  {{0x00741147,0x0151b219,0x01092690,0x00e877e6,0x01f4d6bb,0x0072a332,0x01cd3b03,0x00dadff2,0x0097db5e,0x0086598d},{0x01c69a2b,0x01decf1b,0x02c2fa6e,0x013b7c4f,0x037beac8,0x013a16b5,0x028e7bda,0x01f6e8ac,0x01e34fe9,0x01726947},{0x01f10e67,0x003c73de,0x022b7ea2,0x010f32c2,0x03ff776a,0x00142277,0x01d38b88,0x00776138,0x03c60822,0x01201140}},
	  {{0x0236d175,0x0008748e,0x03c6476d,0x013f4cdc,0x02eed02a,0x00838a47,0x032e7210,0x018bcbb3,0x00858de4,0x01dc7826},{0x00a37fc7,0x0127b40b,0x01957884,0x011d30ad,0x02816683,0x016e0e23,0x00b76be4,0x012db115,0x02516506,0x0154ce62},{0x00451edf,0x00bd749e,0x03997342,0x01cc2c4c,0x00eb6975,0x01a59508,0x03a516cf,0x00c228ef,0x0168ff5a,0x01697b47}},
	  {{0x00527359,0x01783156,0x03afd75c,0x00ce56dc,0x00e4b970,0x001cabe9,0x029e0f6d,0x0188850c,0x0135fefd,0x00066d80},{0x02150e83,0x01448abf,0x02bb0232,0x012bf259,0x033c8268,0x00711e20,0x03fc148f,0x005e0e70,0x017d8bf9,0x0112b2e2},{0x02134b83,0x001a0517,0x0182c3cc,0x00792182,0x0313d799,0x001a3ed7,0x0344547e,0x01f24a0d,0x03de6ad2,0x00543127}},
	  {{0x00dca868,0x00618f27,0x015a1709,0x00ddc38a,0x0320fd13,0x0036168d,0x0371ab06,0x01783fc7,0x0391e05f,0x01e29b5d},{0x01471138,0x00fca542,0x00ca31cf,0x01ca7bad,0x0175bfbc,0x01a708ad,0x03bce212,0x01244215,0x0075bb99,0x01acad68},{0x03a0b976,0x01dc12d1,0x011aab17,0x00aba0ba,0x029806cd,0x0142f590,0x018fd8ea,0x01a01545,0x03c4ad55,0x01c971ff}},
	  {{0x00d098c0,0x000afdc7,0x006cd230,0x01276af3,0x03f905b2,0x0102994c,0x002eb8a4,0x015cfbeb,0x025f855f,0x01335518},{0x01cf99b2,0x0099c574,0x01a69c88,0x00881510,0x01cd4b54,0x0112109f,0x008abdc5,0x0074647a,0x0277cb1f,0x01e53324},{0x02ac5053,0x01b109b0,0x024b095e,0x016997b3,0x02f26bb6,0x00311021,0x00197885,0x01d0a55a,0x03b6fcc8,0x01c020d5}},
	  {{0x02584a34,0x00e7eee0,0x03257a03,0x011e95a3,0x011ead91,0x00536202,0x00b1ce24,0x008516c6,0x03669d6d,0x004ea4a8},{0x00773f01,0x0019c9ce,0x019f6171,0x01d4afde,0x02e33323,0x01ad29b6,0x02ead1dc,0x01ed51a5,0x01851ad0,0x001bbdfa},{0x00577de5,0x00ddc730,0x038b9952,0x00f281ae,0x01d50390,0x0002e071,0x000780ec,0x010d448d,0x01f8a2af,0x00f0a5b7}},
	  {{0x031f2541,0x00d34bae,0x0323ff9d,0x003a056d,0x02e25443,0x00a1ad05,0x00d1bee8,0x002f7f8e,0x03007477,0x002a24b1},{0x0114a713,0x01457e76,0x032255d5,0x01cc647f,0x02a4bdef,0x0153d730,0x00118bcf,0x00f755ff,0x013490c7,0x01ea674e},{0x02bda3e8,0x00bb490d,0x00f291ea,0x000abf40,0x01dea321,0x002f9ce0,0x00b2b193,0x00fa54b5,0x0128302f,0x00a19d8b}},
	  {{0x022ef5bd,0x01638af3,0x038c6f8a,0x01a33a3d,0x039261b2,0x01bb89b8,0x010bcf9d,0x00cf42a9,0x023d6f17,0x01da1bca},{0x00e35b25,0x000d824f,0x0152e9cf,0x00ed935d,0x020b8460,0x01c7b83f,0x00c969e5,0x01a74198,0x0046a9d9,0x00cbc768},{0x01597c6a,0x0144a99b,0x00a57551,0x0018269c,0x023c464c,0x0009b022,0x00ee39e1,0x0114c7f2,0x038a9ad2,0x01584c17}},
	  {{0x03b0c0d5,0x00b30a39,0x038a6ce4,0x01ded83a,0x01c277a6,0x01010a61,0x0346d3eb,0x018d995e,0x02f2c57c,0x000c286b},{0x0092aed1,0x0125e37b,0x027ca201,0x001a6b6b,0x03290f55,0x0047ba48,0x018d916c,0x01a59062,0x013e35d4,0x0002abb1},{0x003ad2aa,0x007ddcc0,0x00c10f76,0x0001590b,0x002cfca6,0x000ed23e,0x00ee4329,0x00900f04,0x01c24065,0x0082fa70}},
	  {{0x02025e60,0x003912b8,0x0327041c,0x017e5ee5,0x02c0ecec,0x015a0d1c,0x02b1ce7c,0x0062220b,0x0145067e,0x01a5d931},{0x009673a6,0x00e1f609,0x00927c2a,0x016faa37,0x01650ef0,0x016f63b5,0x03cd40e1,0x003bc38f,0x0361f0ac,0x01d42acc},{0x02f81037,0x008ca0e8,0x017e23d1,0x011debfe,0x01bcbb68,0x002e2563,0x03e8add6,0x000816e5,0x03fb7075,0x0153e5ac}},
	  {{0x02b11ecd,0x016bf185,0x008f22ef,0x00e7d2bb,0x0225d92e,0x00ece785,0x00508873,0x017e16f5,0x01fbe85d,0x01e39a0e},{0x01669279,0x017c810a,0x024941f5,0x0023ebeb,0x00eb7688,0x005760f1,0x02ca4146,0x0073cde7,0x0052bb75,0x00f5ffa7},{0x03b8856b,0x00cb7dcd,0x02f14e06,0x001820d0,0x01d74175,0x00e59e22,0x03fba550,0x00484641,0x03350088,0x01c3c9a3}},
	  {{0x00dcf355,0x0104481c,0x0022e464,0x01f73fe7,0x00e03325,0x0152b698,0x02ef769a,0x00973663,0x00039b8c,0x0101395b},{0x01805f47,0x019160ec,0x03832cd0,0x008b06eb,0x03d4d717,0x004cb006,0x03a75b8f,0x013b3d30,0x01cfad88,0x01f034d1},{0x0078338a,0x01c7d2e3,0x02bc2b23,0x018b3f05,0x0280d9aa,0x005f3d44,0x0220a95a,0x00eeeb97,0x0362aaec,0x00835d51}},
	  {{0x01b9f543,0x013fac4d,0x02ad93ae,0x018ef464,0x0212cdf7,0x01138ba9,0x011583ab,0x019c3d26,0x028790b4,0x00e2e2b6},{0x033bb758,0x01f0dbf1,0x03734bd1,0x0129b1e5,0x02b3950e,0x003bc922,0x01a53ec8,0x018c5532,0x006f3cee,0x00ae3c79},{0x0351f95d,0x0012a737,0x03d596b8,0x017658fe,0x00ace54a,0x008b66da,0x0036c599,0x012a63a2,0x032ceba1,0x00126bac}},
	  {{0x03dcfe7e,0x019f4f18,0x01c81aee,0x0044bc2b,0x00827165,0x014f7c13,0x03b430f0,0x00bf96cc,0x020c8d62,0x01471997},{0x01fc7931,0x001f42dd,0x00ba754a,0x005bd339,0x003fbe49,0x016b3930,0x012a159c,0x009f83b0,0x03530f67,0x01e57b85},{0x02ecbd81,0x0096c294,0x01fce4a9,0x017701a5,0x0175047d,0x00ee4a31,0x012686e5,0x008efcd4,0x0349dc54,0x01b3466f}},
	  {{0x02179ca3,0x01d86414,0x03f0afd0,0x00305964,0x015c7428,0x0099711e,0x015d5442,0x00c71014,0x01b40b2e,0x01d483cf},{0x01afc386,0x01984859,0x036203ff,0x0045c6a8,0x0020a8aa,0x00990baa,0x03313f10,0x007ceede,0x027429e4,0x017806ce},{0x039357a1,0x0142f8f4,0x0294a7b6,0x00eaccf4,0x0259edb3,0x01311e6e,0x004d326f,0x0130c346,0x01ccef3c,0x01c424b2}},
	  {{0x0364918c,0x00148fc0,0x01638a7b,0x01a1fd5b,0x028ad013,0x0081e5a4,0x01a54f33,0x0174e101,0x003d0257,0x003a856c},{0x00051dcf,0x00f62b1d,0x0143d0ad,0x0042adbd,0x000fda90,0x01743ceb,0x0173e5e4,0x017bc749,0x03b7137a,0x0105ce96},{0x00f9218a,0x015b8c7c,0x00e102f8,0x0158d7e2,0x0169a5b8,0x00b2f176,0x018b347a,0x014cfef2,0x0214a4e3,0x017f1595}},
	  {{0x006d7ae5,0x0195c371,0x0391e26d,0x0062a7c6,0x003f42ab,0x010dad86,0x024f8198,0x01542b2a,0x0014c454,0x0189c471},{0x0390988e,0x00b8799d,0x02e44912,0x0078e2e6,0x00075654,0x01923eed,0x0040cd72,0x00a37c76,0x0009d466,0x00c8531d},{0x02651770,0x00609d01,0x0286c265,0x0134513c,0x00ee9281,0x005d223c,0x035c760c,0x00679b36,0x0073ecb8,0x016faa50}},
	  {{0x02c89be4,0x016fc244,0x02f38c83,0x018beb72,0x02b3ce2c,0x0097b065,0x034f017b,0x01dd957f,0x00148f61,0x00eab357},{0x0343d2f8,0x003398fc,0x011e368e,0x00782a1f,0x00019eea,0x00117b6f,0x0128d0d1,0x01a5e6bb,0x01944f1b,0x012b41e1},{0x03318301,0x018ecd30,0x0104d0b1,0x0038398b,0x03726701,0x019da88c,0x002d9769,0x00a7a681,0x031d9028,0x00ebfc32}},
	  {{0x0220405e,0x0171face,0x02d930f8,0x017f6d6a,0x023b8c47,0x0129d5f9,0x02972456,0x00a3a524,0x006f4cd2,0x004439fa},{0x00c53505,0x0190c2fd,0x00507244,0x009930f9,0x01a39270,0x01d327c6,0x0399bc47,0x01cfe13d,0x0332bd99,0x00b33e7d},{0x0203f5e4,0x003627b5,0x00018af8,0x01478581,0x004a2218,0x002e3bb7,0x039384d0,0x0146ea62,0x020b9693,0x0017155f}},
	  {{0x03c97e6f,0x00738c47,0x03b5db1f,0x01808fcf,0x01e8fc98,0x01ed25dd,0x01bf5045,0x00eb5c2b,0x0178fe98,0x01b85530},{0x01c20eb0,0x01aeec22,0x030b9eee,0x01b7d07e,0x0187e16f,0x014421fb,0x009fa731,0x0040b6d7,0x00841861,0x00a27fbc},{0x02d69abf,0x0058cdbf,0x0129f9ec,0x013c19ae,0x026c5b93,0x013a7fe7,0x004bb2ba,0x0063226f,0x002a95ca,0x01abefd9}},
	  {{0x02f5d2c1,0x00378318,0x03734fb5,0x01258073,0x0263f0f6,0x01ad70e0,0x01b56d06,0x01188fbd,0x011b9503,0x0036d2e1},{0x0113a8cc,0x01541c3e,0x02ac2bbc,0x01d95867,0x01f47459,0x00ead489,0x00ab5b48,0x01db3b45,0x00edb801,0x004b024f},{0x00b8190f,0x011fe4c2,0x00621f82,0x010508d7,0x001a5a76,0x00c7d7fd,0x03aab96d,0x019cd9dc,0x019c6635,0x00ceaa1e}},
	  {{0x01085cf2,0x01fd47af,0x03e3f5e1,0x004b3e99,0x01e3d46a,0x0060033c,0x015ff0a8,0x0150cdd8,0x029e8e21,0x008cf1bc},{0x00156cb1,0x003d623f,0x01a4f069,0x00d8d053,0x01b68aea,0x01ca5ab6,0x0316ae43,0x0134dc44,0x001c8d58,0x0084b343},{0x0318c781,0x0135441f,0x03a51a5e,0x019293f4,0x0048bb37,0x013d3341,0x0143151e,0x019c74e1,0x00911914,0x0076ddde}},
	  {{0x006bc26f,0x00d48e5f,0x00227bbe,0x00629ea8,0x01ea5f8b,0x0179a330,0x027a1d5f,0x01bf8f8e,0x02d26e2a,0x00c6b65e},{0x01701ab6,0x0051da77,0x01b4b667,0x00a0ce7c,0x038ae37b,0x012ac852,0x03a0b0fe,0x0097c2bb,0x00a017d2,0x01eb8b2a},{0x0120b962,0x0005fb42,0x0353b6fd,0x0061f8ce,0x007a1463,0x01560a64,0x00e0a792,0x01907c92,0x013a6622,0x007b47f1}}
	};
#endif
	ge25519 d1;
	ge25519_p1p1 t;
	int32_t i;

	contract256_slidingwindow_modm(slide1, s1, S1_SWINDOWSIZE);
	contract256_slidingwindow_modm(slide2, s2, S2_SWINDOWSIZE);

	ge25519_double(&d1, p1);
	ge25519_full_to_pniels(pre1, p1);
	for (i = 0; i < S1_TABLE_SIZE - 1; i++)
		ge25519_pnielsadd(&pre1[i+1], &d1, &pre1[i]);

	ge25519_set_neutral(r);

	i = 255;
	while ((i >= 0) && !(slide1[i] | slide2[i]))
		i--;

	for (; i >= 0; i--) {
		ge25519_double_p1p1(&t, r);

		if (slide1[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_pnielsadd_p1p1(&t, r, &pre1[abs(slide1[i]) / 2], (unsigned char)slide1[i] >> 7);
		}

		if (slide2[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_nielsadd2_p1p1(&t, r, &table2[abs(slide2[i]) / 2], (unsigned char)slide2[i] >> 7);
		}

		ge25519_p1p1_to_partial(r, &t);
	}
	explicit_bzero(slide1, sizeof(slide1));
	explicit_bzero(slide2, sizeof(slide2));
}

/* computes [s1]p1 */
STATIC void ge25519_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const bignum256modm s1) {
	signed char slide1[256];
	ge25519_pniels pre1[S1_TABLE_SIZE];
	ge25519 d1;
	ge25519_p1p1 t;
	int32_t i;

	contract256_slidingwindow_modm(slide1, s1, S1_SWINDOWSIZE);

	ge25519_double(&d1, p1);
	ge25519_full_to_pniels(pre1, p1);
	for (i = 0; i < S1_TABLE_SIZE - 1; i++)
		ge25519_pnielsadd(&pre1[i+1], &d1, &pre1[i]);

	/* set neutral */
	ge25519_set_neutral(r);

	i = 255;
	while ((i >= 0) && !slide1[i])
		i--;

	for (; i >= 0; i--) {
		ge25519_double_p1p1(&t, r);

		if (slide1[i]) {
			ge25519_p1p1_to_full(r, &t);
			ge25519_pnielsadd_p1p1(&t, r, &pre1[abs(slide1[i]) / 2], (unsigned char)slide1[i] >> 7);
		}

		ge25519_p1p1_to_partial(r, &t);
	}
	explicit_bzero(slide1, sizeof(slide1));
}

/*
 * The following conditional move stuff uses conditional moves.
 * I will check on which compilers this works, and provide suitable
 * workarounds for those where it doesn't.
 *
 * This works on gcc 4.x and above with -O3.  Don't use -O2, this will
 * cause the code to not generate conditional moves.  Don't use any -march=
 * with less than i686 on x86
 */
DONNA_INLINE static void ge25519_cmove_stride4(long * r, long * p, long * pos, long * n, int stride) {
  int i;
  long x0=r[0], x1=r[1], x2=r[2], x3=r[3], y0, y1, y2, y3;
  for(; p<n; p+=stride) {
    int flag=(p==pos);
    y0 = p[0];
    y1 = p[1];
    y2 = p[2];
    y3 = p[3];
    x0 = flag ? y0 : x0;
    x1 = flag ? y1 : x1;
    x2 = flag ? y2 : x2;
    x3 = flag ? y3 : x3;
  }
  r[0] = x0;
  r[1] = x1;
  r[2] = x2;
  r[3] = x3;
}
#define HAS_CMOVE_STRIDE4

DONNA_INLINE static void ge25519_cmove_stride4b(long * r, long * p, long * pos, long * n, int stride) {
  int i;
  long x0=p[0], x1=p[1], x2=p[2], x3=p[3], y0, y1, y2, y3;
  for(p+=stride; p<n; p+=stride) {
    int flag=(p==pos);
    y0 = p[0];
    y1 = p[1];
    y2 = p[2];
    y3 = p[3];
    x0 = flag ? y0 : x0;
    x1 = flag ? y1 : x1;
    x2 = flag ? y2 : x2;
    x3 = flag ? y3 : x3;
  }
  r[0] = x0;
  r[1] = x1;
  r[2] = x2;
  r[3] = x3;
}
#define HAS_CMOVE_STRIDE4B

STATIC void ge25519_move_conditional_pniels_array(ge25519_pniels * r, const ge25519_pniels * p, int pos, int n) {
#ifdef HAS_CMOVE_STRIDE4B
  int i;
  for(i=0; i<sizeof(ge25519_pniels)/sizeof(long); i+=4) {
    ge25519_cmove_stride4b(((long*)r)+i,
			   ((long*)p)+i,
			   ((long*)(p+pos))+i,
			   ((long*)(p+n))+i,
			   sizeof(ge25519_pniels)/sizeof(long));
  }
#else
  int i;
  for(i=0; i<n; i++) {
    ge25519_move_conditional_pniels(r, p+i, pos==i);
  }
#endif
}

STATIC void ge25519_move_conditional_niels_array(ge25519_niels * r, const uint8_t p[8][96], int pos, int n) {
  int i;
  for(i=0; i<96/sizeof(long); i+=4) {
    ge25519_cmove_stride4(((long*)r)+i,
			  ((long*)p)+i,
			  ((long*)(p+pos))+i,
			  ((long*)(p+n))+i,
			  96/sizeof(long));
  }
}

/* computes [s1]p1, constant time */
STATIC void ge25519_scalarmult(ge25519 *r, const ge25519 *p1, const bignum256modm s1) {
	signed char slide1[64];
	ge25519_pniels pre1[9];
	ge25519_pniels pre;
	ge25519 d1, r1;
	ge25519_p1p1 t;
	int32_t i, j;

	contract256_window4_modm(slide1, s1);

	/* set neutral */
	ge25519_set_neutral(r);

	ge25519_full_to_pniels(pre1, r);
	ge25519_full_to_pniels(pre1+1, p1);
	ge25519_double(&d1, p1);
	ge25519_full_to_pniels(pre1+2, &d1);
	for (i = 1; i < 7; i++) {
		ge25519_pnielsadd(&pre1[i+2], &d1, &pre1[i]);
	}

	for (i = 63; i >= 0; i--) {
		int k=abs(slide1[i]);
		ge25519_double_partial(r, r);
		ge25519_double_partial(r, r);
		ge25519_double_partial(r, r);
		ge25519_double_p1p1(&t, r);
		ge25519_move_conditional_pniels_array(&pre, pre1, k, 9);
		ge25519_p1p1_to_full(r, &t);
		ge25519_pnielsadd_p1p1(&t, r, &pre, (unsigned char)slide1[i] >> 7);
		ge25519_p1p1_to_partial(r, &t);
	}
	explicit_bzero(slide1, sizeof(slide1));
}

#if !defined(HAVE_GE25519_SCALARMULT_BASE_CHOOSE_NIELS)

DONNA_INLINE static uint32_t
ge25519_windowb_equal(uint32_t b, uint32_t c) {
	return ((b ^ c) - 1) >> 31;
}

#include <stdio.h>
static void
ge25519_scalarmult_base_choose_niels(ge25519_niels *t, const uint8_t table[256][96], uint32_t pos, signed char b) {
	bignum25519 neg;
	uint32_t sign = (uint32_t)((unsigned char)b >> 7);
	uint32_t mask = ~(sign - 1);
	uint32_t u = (b + mask) ^ mask;
	uint32_t i;

	/* ysubx, xaddy, t2d in packed form. initialize to ysubx = 1, xaddy = 1, t2d = 0 */
	uint8_t packed[96] = {0};
	packed[0] = 1;
	packed[32] = 1;

	ge25519_move_conditional_niels_array(packed, &table[pos*8], u-1, 8);

	/* expand in to t */
	curve25519_expand(t->ysubx, packed +  0);
	curve25519_expand(t->xaddy, packed + 32);
	curve25519_expand(t->t2d  , packed + 64);

	/* adjust for sign */
	curve25519_swap_conditional(t->ysubx, t->xaddy, sign);	
	curve25519_neg(neg, t->t2d);
	curve25519_swap_conditional(t->t2d, neg, sign);
}

#endif /* HAVE_GE25519_SCALARMULT_BASE_CHOOSE_NIELS */


/* computes [s]basepoint */
static void
ge25519_scalarmult_base_niels(ge25519 *r, const uint8_t basepoint_table[256][96], const bignum256modm s) {
	signed char b[64];
	uint32_t i;
	ge25519_niels t;

	contract256_window4_modm(b, s);

	ge25519_scalarmult_base_choose_niels(&t, basepoint_table, 0, b[1]);
	curve25519_sub_reduce(r->x, t.xaddy, t.ysubx);
	curve25519_add_reduce(r->y, t.xaddy, t.ysubx);
	memset(r->z, 0, sizeof(bignum25519));
	curve25519_copy(r->t, t.t2d);
	r->z[0] = 2;	
	for (i = 3; i < 64; i += 2) {
		ge25519_scalarmult_base_choose_niels(&t, basepoint_table, i / 2, b[i]);
		ge25519_nielsadd2(r, &t);
	}
	ge25519_double_partial(r, r);
	ge25519_double_partial(r, r);
	ge25519_double_partial(r, r);
	ge25519_double(r, r);
	ge25519_scalarmult_base_choose_niels(&t, basepoint_table, 0, b[0]);
	curve25519_mul(t.t2d, t.t2d, ge25519_ecd);
	ge25519_nielsadd2(r, &t);
	for(i = 2; i < 64; i += 2) {
		ge25519_scalarmult_base_choose_niels(&t, basepoint_table, i / 2, b[i]);
		ge25519_nielsadd2(r, &t);
	}
}

STATIC void ge25519_scalarmult_base(ge25519 *r, const bignum256modm s) {
	ge25519_scalarmult_base_niels(r, ge25519_niels_base_multiples, s);
}
