#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <gmp.h>
#include <sys/time.h>
#include <string.h>
#include "sha256.c"
#include <ctype.h>

#define INF 1
#define nPoints 20

struct coord
{
	mpz_t x;
	mpz_t y;
	int inf;
};

struct parameters{
	mpz_t p;
 	mpz_t a;
	mpz_t b;
	struct coord G;
	mpz_t n;
	// unsigned int h; // não é útil
};

struct pe
{
	struct coord pu;
	struct coord id;
	struct coord ids;
	struct coord V;
	struct coord R;
	mpz_t k;
	mpz_t k_sess;
	mpz_t v;
	mpz_t r;
	mpz_t sa;
	mpz_t s;
};

int mSize;
struct parameters ec; // parâmetros globais
struct coord points[nPoints]; //pontos[índice do char]
struct pe alice; // chave pública e privada
struct pe bob; // chave pública e privada
mpz_t seed; // random seed
gmp_randstate_t state; // random state

double timestamp(void);

void showPoints();
void showPoint(struct coord p);

void modn(mpz_t rop, mpz_t a);
void modp(mpz_t rop, mpz_t a);

void eccDbl(struct coord *rop, struct coord p);
void eccAdd(struct coord *rop, struct coord p, struct coord q);
void eccSub(struct coord *rop, struct coord p, struct coord q);
void mult(struct coord *rop, mpz_t k, struct coord p);

void euclidian(mpz_t rop, mpz_t a, mpz_t md);
void multInv(mpz_t rop, mpz_t a, mpz_t md);

void generateKeys();
bool isValidPoint(mpz_t x, mpz_t y);
void findPoints();

void eccCipher(struct coord *c2, struct coord *c1, int size);
void eccDecipher(struct coord *d, struct coord *c2, struct coord c1, int s);

void test(long int messageSize);

void print_hash(unsigned char hash[]);
void hash(mpz_t rop, mpz_t m);

void generateCertificates();
void sign(mpz_t m);
int check(mpz_t m);