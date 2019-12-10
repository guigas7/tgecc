#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "sha256.c"


//~ Tamanho do P. Também ditará a ordem da curva.
#define P 23 

// Parâmetro A da curva.
//~ #define A 43
#define A 1

//~ Parâmetro B da curva.
//~ #define B 5
#define B 4

//~ Número de PEs no sistema local
#define N 3

//~ Personagens para teste
#define PKG		0
#define ALICE	1
#define BOB		2

/* Ordem da Curva Elíptica.
 * O Teorema de Hasse estabelece que o número de pontos na curva é chamado de ordem da curva elíptica.
 * Esta ordem é roughly igual ao tamanho de p (Zp) do corpo finito.
 */
int nPoints = -1;


struct coord
{
	int x;
	int y;
};
struct pe
{
	struct coord pu;
	int k;
};

struct coord *points = NULL;
struct coord G;
struct pe PEs[N];
char alpha[26] = "abcdefghijklmnopqrstuvxwyz";


//~ Operação modular
int mod(int a)
{
	return (a < 0) ? (a += P) % P : a % P;
}

//~ Calcular o número de pontos na curva
int numberPoints()
{
	//~ int min = ceil(P + 1 - 2 * sqrt(P));
	int max = floor(P + 1 + 2 * sqrt(P));
	//~ printf("Number of points = %.0f < Np < %.0f\n",
		   //~ ceil((P + 1) - 2 * sqrt(P)), floor(P + 1 + 2 * sqrt(P)));
	return max;
}

//~ Algoritmo extendido de euclides
int extEuclid(int a)
{
	int u = 1;
	int v = 0;
	int g = a;
	int u1 = 0;
	int v1 = 1;
	int g1 = P;
	while (g1 != 0)
	{
		int q = g / g1;
		int t1 = u - q * u1;
		int t2 = v - q * v1;
		int t3 = g - q * g1;
		u = u1;
		v = v1;
		g = g1;
		u1 = t1;
		v1 = t2;
		g1 = t3;
	}
	return mod(u);
}

//~ Inveso Multiplicativo
int multInv(int a)
{
	return extEuclid(a);
}

//~ ECCADD r=p+q
struct coord eccAdd(struct coord p, struct coord q)
{
	struct coord r;
	int s = mod(mod(q.y - p.y) * mod(multInv(mod(q.x - p.x))));	// Slope
	r.x = mod(mod(mod(s * s) - p.x) - q.x);
	r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
	r.x = mod(r.x);
	r.y = mod(r.y);
	return r;
}

//~  r=2p
struct coord eccDbl(struct coord p)
{
	struct coord r;
	if (p.y == -p.y || p.y == 0) {
		r.x = r.y = 0;
	}
	else if (p.y != 0)
	{
		int s = mod(mod(3 * (p.x * p.x) + A) * mod(multInv(mod(2 * p.y))));
		printf("mod(s * s) - mod(p.x * 2): %d\n", mod(s * s) - mod(p.x * 2));
		printf("mod(mod(s * s) - mod(p.x * 2)): %d\n", mod(mod(s * s) - mod(p.x * 2)));
		r.x = mod(mod(s * s) - mod(2 * p.x));
		r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
	}
	return r;
}

//~ ECCSUB r=p-q
struct coord eccSub(struct coord p, struct coord q)
{
	struct coord q1;
	q1.x = q.x;
	q1.y = mod(q.y * -1);
	struct coord r = eccAdd(p, q1);
	return r;
}

//~ Left-to-right binary algorithm
struct coord mult(int k, struct coord p)
{
	//~ output: q = kP
	struct coord r0 = p;
	struct coord r1 = p;
	int n = k;
	for (int i = n - 2; i > 0; i--)
	{
		if (i > 1)
			r0 = eccDbl(r0);
		else if (i == 1)
			r0 = eccAdd(r1, r0);
	}
	return r0;
}

//~ Verifica a validade do ponto
bool isValidPoint(int x, int y)
{
	if (mod(y * y) == mod(mod(x * x * x) + mod(A * x) + B))
		return true;
	return false;
}
void showPoints()
{
	printf("\nPoints:\n");
	for (int i = 0; i < nPoints; i++)
		printf("%.2d,%.2d\t", points[i].x, points[i].y);
	printf("\n");
}

void showPoint(struct coord p)
{
	printf("\nPoint:\n");
	printf("%.2d,%.2d", p.x, p.y);
	printf("\n");
}
//~ Cálculo de todos pontos da curva
void calculatePoints()
{
	nPoints = numberPoints();
	points = malloc(sizeof(struct coord) * nPoints);
	int k = 0;
	for (int i = 0; i < P; i++)
		for (int j = 0; j < P; j++)
			if (isValidPoint(i, j))
			{
				points[k].x = i;
				points[k].y = j;
				k++;
			}
	nPoints = k;
	//~ printf("%d generated points + Infinity Point (0,1,0)\n", k);
}

//~ Selecione aleatoriamente um ponto gerador (G)
void generateG()
{
	srand(time(NULL));			//~ Mudar a semente   
	G = points[P - (rand() % P)];
	showPoint(G);
}

//~ Gerar as chaves públicas e privadas
void generateKeys()
{
	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		PEs[i].k = 1 + (rand() % P);	//~ Chave privada (k)
		PEs[i].pu = mult(PEs[i].k, G);	//~ Chaves Públicas
	}
}

void generateMessage(int size, struct coord *m, char *plaintext)
{
	for (int i = 0; i < size; i++)
	{
		int index = (rand() % 25);
		plaintext[i] = alpha[index];
		m[i].x = points[index].x;
		m[i].y = points[index].y;
	}
	//~ printf("\nGenerated message: %s\n", plaintext);
	printf("\nmessage\n");
	for (int i = 0; i < size; i++)
		printf("(%d,%d)\n", m[i].x, m[i].y);

}

void eccCipher(struct pe a, struct pe b, struct coord *m, struct coord *c,
			   int s)
{
	for (int i = 0; i < s; i++)
		c[i] = eccAdd(m[i], mult(a.k, b.pu));

	printf("\n\nGenerated ciphered message in points: \n");
	for (int i = 0; i < s; i++)
		printf("(%d,%d)\n", c[i].x, c[i].y);
}

void eccDecipher(struct pe a, struct pe b, struct coord *d, struct coord *c,
				 int s)
{
	for (int i = 0; i < s; i++)
		d[i] = eccSub(c[i], mult(b.k, a.pu));

	printf("\ndeciphered\n");
	printf("\n\nGenerated deciphered message in points: \n");
	for (int i = 0; i < s; i++)
		printf("(%d,%d)\n", d[i].x, d[i].y);
	printf("\n");
}

int hash(char *m, int size)
{
	printf("\n - m:\t\t%s (%d bytes)\n", m, size);

	//~ SHA256 convertido para inteiro modular
	int result = 0;
	char *hash = malloc(sizeof(char) * size);
	SHA256_CTX ctx;
	sha256_init(&ctx);
	sha256_update(&ctx, m, strlen(m));
	sha256_final(&ctx, hash, &result, size, P);

	result = mod(result);
	printf("\n - Int hash:\t%d \n", result);

	return result;
}

//~ ECDSA SIGN
struct coord sign(int d, char *m, int size)
{
	int h = hash(m, size);
	srand(time(NULL));
	int r = 0, s = 0;
	do
	{
		int k = 1 + (rand() % P);
		struct coord kG = mult(k, G);
		r = mod(kG.x);
		s = mod(multInv(k) * (h + r * d));
	}
	while (r == 0 || s == 0);

	struct coord sign;
	sign.x = r;
	sign.y = s;
	printf("\nGenerated signature in points: (%d,%d) ", r, s);
	return sign;
}

bool checkSignature(struct coord Q, struct coord sign, char *m)
{
	int r = sign.x;
	int s = sign.y;

	printf("\nSignature Checking:");
	printf("\n - r:\t\t%d", r);
	printf("\n - s:\t\t%d", s);

	if (r >= 1 && s >= 1 && r <= (P - 1) && s <= (P - 1))
	{
		int w = mod(multInv(s));
		int h = hash(m, strlen(m));
		int u1 = mod(h * w);
		int u2 = mod(r * w);

		printf(" - w:\t\t%d", w);
		printf("\n - h:\t\t%d (int hash)", h);
		printf("\n - u1:\t\t%d", u1);
		printf("\n - u2:\t\t%d", u2);

		struct coord u1G = mult(u1, G);
		struct coord u2Q = mult(u2, Q);
		struct coord v = eccAdd(u1G, u2Q);

		printf("\n - u1G:\t\t%d(%d,%d) = (%d,%d)", u1, G.x, G.y, u1G.x, u1G.y);
		printf("\n - u2Q:\t\t%d(%d,%d) = (%d,%d)", u2, Q.x, Q.y, u2Q.x, u2Q.y);
		printf("\n - v: = u1G+u2Q:\t(%d,%d) + (%d,%d) = (%d,%d)", u1G.x, u1G.y,
			   u2Q.x, u2Q.y, v.x, v.y);
		printf("\n - Result (r == v.x): \t(%d == %d): ", r, v.x);
		if (r == v.x)
			return true;
	}
	return false;
}

//~ It's show time
void main()
{
	//double x = pow(2,192) - pow(2,32) - pow(2,12) - pow(2,8) - pow(2,7) - pow(2,6) - pow(2,3) - 1;
	//printf("%lf",x);

	// calculatePoints();
	// generateG();
	// showPoints();
	// generateKeys();
	// printf("\n - G:\t\t[%d,%d]\n", G.x, G.y);
	// printf(" - PU[kgc]:\t[%d,%d]\n", PEs[PKG].pu.x, PEs[PKG].pu.y);
	// printf(" - PU[Alice]:\t[%d,%d]\n", PEs[ALICE].pu.x, PEs[ALICE].pu.y);
	// printf(" - PU[Bob]:\t[%d,%d]\n", PEs[BOB].pu.x, PEs[BOB].pu.y);
	G.x = 1;
	G.y = 12;
	struct coord teste = eccDbl(G);
	printf("teste.x: %d, teste.y: %d\n", teste.x, teste.y);
	if (isValidPoint(teste.x, teste.y)) {
		printf("deu caralho22!\n");
	}

	// // Mensagem
	// int messageSize = 32;		// bytes
	// char *plaintext = malloc(sizeof(char) * messageSize);

	// struct coord message[messageSize];
	// generateMessage(messageSize, message, plaintext);

	// // Cifra
	// struct coord ciphered[messageSize];
	// eccCipher(PEs[ALICE], PEs[BOB], message, ciphered, messageSize);

	// // Decifra
	// struct coord deciphered[messageSize];
	// eccDecipher(PEs[ALICE], PEs[BOB], deciphered, ciphered, messageSize);

	// // Assina  
	// struct coord signature = sign(PEs[ALICE].k, plaintext, messageSize);

	// // Verifica a assinatura
	// if (checkSignature(PEs[ALICE].pu, signature, plaintext))
	// 	printf("Assinatura válida!");
	// else
	// 	printf("Assinatura inválida!");

	// printf("\n\n");
	// points = NULL;
}
