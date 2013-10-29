#ifndef ED25519_PRIMS_H
#define ED25519_PRIMS_H

void ge25519_pack(unsigned char r[32], const ge25519 *p);
int ge25519_unpack_negative_vartime(ge25519 *r, const unsigned char p[32]);
void ge25519_double_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const bignum256modm s1, const bignum256modm s2);
void ge25519_scalarmult_vartime(ge25519 *r, const ge25519 *p1, const bignum256modm s1);
void ge25519_scalarmult(ge25519 *r, const ge25519 *p1, const bignum256modm s1);
void ge25519_scalarmult_base(ge25519 *r, const bignum256modm s);
void expand256_modm(bignum256modm out, const unsigned char *in, size_t len);
void expand_raw256_modm(bignum256modm out, const unsigned char in[32]);
void contract256_modm(unsigned char out[32], const bignum256modm in);
void add256_modm(bignum256modm r, const bignum256modm x, const bignum256modm y);
void mul256_modm(bignum256modm r, const bignum256modm x, const bignum256modm y);
void sub256_modm_batch(bignum256modm out, const bignum256modm a, const bignum256modm b, size_t limbsize);
int lt256_modm_batch(const bignum256modm a, const bignum256modm b, size_t limbsize);
int lte256_modm_batch(const bignum256modm a, const bignum256modm b, size_t limbsize);
int iszero256_modm_batch(const bignum256modm a);
int isone256_modm_batch(const bignum256modm a);

#endif
