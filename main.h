#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include "sha256.c"
#include <omp.h>

//~ Número de PEs no sistema local
#define N 2

//~ Personagens para teste
#define ALICE	0
#define BOB		1

#define pamnt 30

/* Ordem da Curva Elíptica.
 * O Teorema de Hasse estabelece que o número de pontos na curva é chamado de ordem da curva elíptica.
 * Esta ordem é roughly igual ao tamanho de p (Zp) do corpo finito.
 */
mpz_t nPoints;

struct coord
{
	mpz_t x;
	mpz_t y;
};

struct parameters{
	mpz_t p;
 	mpz_t a;
	mpz_t b;
	struct coord G;
	mpz_t n;
	unsigned int h;
};

struct parameters ec;

struct pe
{
	struct coord pu;
	mpz_t k;
};

struct coord points[pamnt];
struct pe PEs[N];
char alpha[26] = "abcdefghijklmnopqrstuvxwyz";

void showPoints();
void showPoint(struct coord p);

void modn(mpz_t res, mpz_t a);
void modp(mpz_t res, mpz_t a);

void eccAdd(struct coord *rop, struct coord p, struct coord q);
void eccDbl(struct coord *rop, struct coord p);
void eccSub(struct coord *rop, struct coord p, struct coord q);
void mult(struct coord *rop, mpz_t k, struct coord p);
void euclidian(mpz_t rop, mpz_t a);
void multInv(mpz_t rop, mpz_t a);

void generateKeys();
void generateMessage(int size, struct coord *m, char *plaintext);

void eccCipher(struct pe a, struct pe b, struct coord *m, struct coord *c, int s);
void eccDecipher(struct pe a, struct pe b, struct coord *d, struct coord *c, int s);

bool isValidPoint(mpz_t x, mpz_t y);
bool isnotduplicate(struct coord c, int k);
void findPoints();

// 	g.x 48439561293906451759052585252797914202762949526041747995844080717082404635286
//	g.y 36134250956749795798585127919587881956611106672985015071877198253568414405109
//	x1. 56515219790691171413109057904011688695424810155802929973526481321309856242040
//	y1. 03377031843712258259223711451491452598088675519751548567112458094635497583569

/*
// não é o do ivan
void calculatePoints()
{
	int size = 1000000;
	pamnt = 30;
	struct coord *vals = malloc(sizeof(struct coord) * size);
	points = malloc(sizeof(struct coord) * pamnt);
	mpz_t aux1, aux2, aux3, nsq;
	mpz_init(aux1);
	mpz_init(aux2);
	mpz_init(aux3);
	mpz_init(nsq);
	mpz_init(n);
	mpz_init(size);
	// usar n como um mpz_t e fazer o controle de alocação
    for (mpz_set_ui(n, 0); n < size; ++n) {
    	printf("%d\n", n);
        // assignment 1 -> points[n].x = ((n*n*n) + ec.a * n + ec.b) % P;
		// nsq = n * n;    	
    	mpz_set_ui(nsq, n * n);
    	// aux1 = n * n * n;
    	mpz_mul_ui(aux1, nsq, n);
    	// aux2 = ec.a * n;
    	mpz_mul_ui(aux2, ec.a, n);
    	// aux3 = (n*n*n) + ec.a * n
    	mpz_add(aux3, aux1, aux2);
    	// aux1 = (n*n*n) + ec.a * n + ec.b
    	mpz_add(aux1, aux3, ec.b);
    	// points[n].x = ((n*n*n) + ec.a * n + ec.b) % P
    	//printf("printou esse\n");
    	mpz_init(vals[n].x);
    	//printf("mas não esse");
    	modp(vals[n].x, aux1);
        // assignment 2 -> points[n].y = (n * n) % P;
    	mpz_init(vals[n].y);
    	modp(vals[n].y, nsq);
    }
	mpz_clear(aux1);
	mpz_clear(aux2);
	mpz_clear(aux3);
	mpz_clear(nsq);
    int k = 0;
    for (int n = 0; n < size; ++n) {
        for (int m = 0; m < size; ++m) {
            if (mpz_cmp(vals[n].x, vals[m].y) == 0) {
            	mpz_init_set_ui(points[k].x, n);
            	mpz_init_set_ui(points[k].y, m);
            	++k;
            	printf("point %d:\nx: %d\ny: %d", k, n, m);
            }
        }
    }
    pamnt = k;
    // for (int n = 0; n < size; ++n) {
    // 	mpz_clears(points[n].x, points[n].y);
    // }
}

//~ Cálculo de todos pontos da curva do ivan
void calculatePoints()
{
	mpz_t i, j;
	mpz_init(nPoints);
	numberPoints(nPoints);
	
	int k = 0;
	// fazer o break para testar a função isvalidpoint
	// #pragma omp parallel for shared(k) collapse(2) private(j) num_threads(3)
	for ( mpz_init_set_ui(i, 0); mpz_cmp(i, ec.p) < 0 && k < pamnt; mpz_add_ui(i, i, 1)) {
		for ( mpz_init_set_ui(j, 0); mpz_cmp(j, ec.p) < 0 && k < pamnt; mpz_add_ui(j, j, 1)) {
			if (isValidPoint(i, j))
			{
				mpz_init_set(points[k].x, i);
				mpz_init_set(points[k].y, j);
				k++;
				gmp_printf("found (i: %Zd, j: %Zd)\n", i, j);
			}
		}
	}
	printf("k: %d\n", k);
	// printf("%d generated points + Infinity Point (0,1,0)\n", k);
}

//~ Cálculo de todos pontos da curva
void calculatePoints()
{
	mpz_t i, j;
	mpz_init(nPoints);
	numberPoints(nPoints);
	pamnt = 30;
	points = malloc(sizeof(struct coord) * pamnt);
	int k = 0;
	// fazer o break para testar a função isvalidpoint
	// #pragma omp parallel for shared(k) collapse(2) private(j) num_threads(3)
	for ( mpz_init_set_ui(i, 0); mpz_cmp(i, ec.p) < 0 && k < pamnt; mpz_add_ui(i, i, 1)) {
		for ( mpz_init_set_ui(j, 0); mpz_cmp(j, ec.p) < 0 && k < pamnt; mpz_add_ui(j, j, 1)) {
			if (isValidPoint(i, j))
			{
				mpz_init_set(points[k].x, i);
				mpz_init_set(points[k].y, j);
				k++;
				gmp_printf("found (i: %Zd, j: %Zd)\n", i, j);
			}
		}
	}
	printf("k: %d\n", k);
	// printf("%d generated points + Infinity Point (0,1,0)\n", k);
}

//~ Calcular o número de pontos na curva
void numberPoints(mpz_t res)
{
	mpf_t fp, aux1, aux2;;
	mpf_set_default_prec(512);
  	mpf_init_set_str(fp, "115792089210356248762697446949407573530086143415290314195533631308867097853951.0", 10);
	mpf_init(aux1);
	mpf_init(aux2);
	// npoints = floor(p + 1 + 2 * sqrt(p));

	//gmp_printf("fp                        : %Ff\n", fp);
	// aux1 = sqrt(p)
	mpf_sqrt(aux1, fp);
	//gmp_printf("sqrt(p)                   : %Ff\n", aux1);
	// aux2 = 2 * sqrt(p)
	mpf_mul_ui(aux2, aux1, 2);
	//gmp_printf("2 * sqrt(p)               : %Ff\n", aux2);
	// aux1 = 1 + 2 * sqrt(p)
	mpf_add_ui(aux1, aux2, 1);
	//gmp_printf("1 + 2 * sqrt(p)           : %Ff\n", aux1);
	// aux2 = p + 1 + 2 * sqrt(p)
	mpf_add(aux2, fp, aux1);
	//gmp_printf("p + 1 + 2 * sqrt(p)       : %Ff\n", aux2);
	// aux1 = floor(p + 1 + 2 * sqrt(p))
	mpf_floor(aux1, aux2);
	//gmp_printf("floor(p + 1 + 2 * sqrt(p)): %Ff\n", aux1);
	mpz_set_f(res, aux1); // returning result
	//gmp_printf("npoints: %Zd\n", res);
	mpf_clear(aux1);
	mpf_clear(aux2);
	mpf_clear(fp);
}

//~ Verifica a validade do ponto
bool isValidPoint(mpz_t x, mpz_t y)
{
	// mod(y * y) == mod(mod(x * x * x) + mod(A * x) + B)
	mpz_t op, aux1, aux2, res;
	mpz_init(aux1);
	mpz_init(aux2);
	mpz_init(op);
	mpz_init(res);
	// gmp_printf("checking: (x: %Zd, y: %Zd)\n", x, y);
	// op = mod(y * y)
	mpz_set(aux1, y);
	mpz_mul(res, aux1, y);
	modp(op, res);
	// gmp_printf("mod(y * y): %Zd\n", op);
	// aux1 = mod(x * x * x)
	mpz_mul(aux1, x, x);
	mpz_mul(res, aux1, x);
	modp(aux1, res);
	// gmp_printf("mod(x * x * x): %Zd\n", aux1);
	// aux2 = mod(A * X)
	mpz_mul(res, ec.a, x);
	modp(aux2, res);
	// gmp_printf("mod(A * x): %Zd\n", op);
	// aux2 = mod(x * x * x) + mod(A * x) + B
	mpz_add(res, aux2, ec.b);
	mpz_add(aux2, aux1, res);
	if (mpz_cmp(op, aux2) == 0) {
		mpz_clear(aux1);
		mpz_clear(aux2);
		mpz_clear(op);
		mpz_clear(res);
		return true;
	}
	mpz_clear(aux1);
	mpz_clear(aux2);
	mpz_clear(op);
	mpz_clear(res);
	return false;
}
*/
