#ifndef ED25519_PRIMS_H
#define ED25519_PRIMS_H

#include <stdlib.h>

/* cpu */
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__ ) || defined(_M_X64)
        #define CPU_X86_64
#elif defined(__i586__) || defined(__i686__) || (defined(_M_IX86) && (_M_IX86 >= 500))
        #define CPU_X86 500
#elif defined(__i486__) || (defined(_M_IX86) && (_M_IX86 >= 400))
        #define CPU_X86 400
#elif defined(__i386__) || (defined(_M_IX86) && (_M_IX86 >= 300)) || defined(__X86__) || defined(_X86_) || defined(__I86__)
        #define CPU_X86 300
#elif defined(__ia64__) || defined(_IA64) || defined(__IA64__) || defined(_M_IA64) || defined(__ia64)
        #define CPU_IA64
#endif

#if defined(__sparc__) || defined(__sparc) || defined(__sparcv9)
        #define CPU_SPARC
        #if defined(__sparcv9)
                #define CPU_SPARC64
        #endif
#endif

#if defined(powerpc) || defined(__PPC__) || defined(__ppc__) || defined(_ARCH_PPC) || defined(__powerpc__) || defined(__powerpc) || defined(POWERPC) || defined(_M_PPC)
        #define CPU_PPC
        #if defined(_ARCH_PWR7)
                #define CPU_POWER7
        #elif defined(__64BIT__)
                #define CPU_PPC64
        #else
                #define CPU_PPC32
        #endif
#endif

#if defined(__hppa__) || defined(__hppa)
        #define CPU_HPPA
#endif

#if defined(__alpha__) || defined(__alpha) || defined(_M_ALPHA)
        #define CPU_ALPHA
#endif

#if defined(CPU_X86_64) || defined(CPU_IA64) || defined(CPU_SPARC64) || defined(__64BIT__) || defined(__LP64__) || defined(_LP64) || (defined(_MIPS_SZLONG) && (_MIPS_SZLONG == 64))
        #define CPU_64BITS
#endif

#if defined(__GNUC__)
        #define COMPILER_GCC
#endif

/* uint128_t */
#if defined(CPU_64BITS)
        #if defined(COMPILER_GCC)
                #define HAVE_NATIVE_UINT128
                #define HAVE_UINT128
                typedef unsigned uint128_t __attribute__((mode(TI)));

                #define mul64x64_128(out,a,b) out = (uint128_t)a * b;
                #define shr128_pair(out,hi,lo,shift) out = (uint64_t)((((uint128_t)hi << 64) | lo) >> shift);
                #define shr128(out,in,shift) out = (uint64_t)(in >> shift);
                #define add128(a,b) a += b;
                #define add128_64(a,b) a += b;
                #define lo128(a) ((uint64_t)a)
        #else
                #warning need 128bit define for this compiler
        #endif
#endif

#if defined(ED25519_SSE2)
typedef uint32_t bignum25519[12];
typedef uint64_t bignum256modm_element_t;
typedef bignum256modm_element_t bignum256modm[5];
#elif defined(HAVE_UINT128)
typedef uint64_t bignum256modm_element_t;
typedef bignum256modm_element_t bignum256modm[5];
typedef uint64_t bignum25519[5];
#else
typedef uint32_t bignum256modm_element_t;
typedef bignum256modm_element_t bignum256modm[9];
typedef uint32_t bignum25519[10];
#endif

typedef struct ge25519_t {
        bignum25519 x, y, z, t;
} ge25519;

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
