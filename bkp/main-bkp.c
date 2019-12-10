#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

// Tamanho do P. Também ditará a ordem da curva.
#define P 23

// Parâmetro A da curva.
#define A 1

// Parâmetro B da curva.
#define B 0

// Número de PEs no sistema local
#define N 3

//~ Personagens para teste
#define PKG		0
#define ALICE	1
#define BOB		2

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
struct coord points[P];
struct coord G;
struct pe PEs[N];
//~ Operação modular
int mod(int a)
{
	int b = (a < 0) ? (a += P) % P : a % P;	
	printf("\nMod: %d mod %d = %d",a,P,b);	
	return b;
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
	return mod(extEuclid(a));
}

//~ ECCADD r=p+q
struct coord eccAdd(struct coord p, struct coord q)
{
	struct coord r;

	// Propriedades do grupo multiplicativo.
	// identidade aditivia
	if (p.x == 0 && p.y == 0)
	{
		r.x = q.x;
		r.y = q.y;
		return r;
	}
	if (q.x == 0 && q.y == 0)
	{
		r.x = p.x;
		r.y = p.y;
		return r;
	}
	if (p.y == -q.y)
	{
		r.x = r.y = 0;
		return r;
	}

	// Slope
	int s = mod(mod(p.y - q.y) * mod(multInv(mod(p.x - q.y))));
	// X
	r.x = mod(mod(mod(s * s) - p.x) - q.x);
	// Y
	if (s != 0)
		r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
	else
		r.x = r.y = 0;

	//~ printf("\n[%.2d,%.2d]+[%.2d,%.2d] -> [%.2d,%.2d]", p.x, p.y, q.x, q.y, r.x, r.y);

	return r;
}

//~ ECCDBL r=2p
struct coord eccDbl(struct coord p)
{
	struct coord r;
	// Propriedades do grupo multiplicativo. 
	// identidade aditivia
	if (p.x == 0 && p.y == 0)
	{
		r.x = p.x;
		r.y = p.y;
		return r;
	}
	if (p.x == 0 && p.y == 0)
	{
		r.x = p.x;
		r.y = p.y;
		return r;
	}
	if (p.y == -p.y)
	{
		r.x = r.y = 0;
		return r;
	}

	// Slope
	int s = (mod(mod(3 * p.x) + A)) * multInv(mod(2 * p.y));
	// X
	r.x = mod(mod(s * s) - mod(2 * p.x));
	// Y
	if (s != 0)
		r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
	else
		r.x = r.y = 0;

	//~ printf("\n[%.2d,%.2d]+[%.2d,%.2d] -> [%.2d,%.2d]", p.x, p.y, p.x, p.y, r.x, r.y);

	return r;
}

//~ Multiplicalção Escalar Binário para Corpos Finitos
struct coord eccMultBinary(int k, struct coord p)
{
	//~ q = kp
	struct coord q = p;
	for (int i = N - 1; i > 0; i--)
	{
		q = eccDbl(q);
		if (i == 1)
			q = eccAdd(q, p);
	}
	return q;
}

//~ Verifica a validade do ponto
bool isValidPoint(int x, int y)
{
	if (mod(y * y) == mod(mod(x * x * x) + mod(A * x) + B))
		return true;
	return false;
}

//~ Cálculo de todos pontos da curva
void calculatePoints()
{
	int k = 0;
	for (int i = 1; i <= P; ++i)
	{
		for (int j = 1; j <= P; ++j)
		{
			if (isValidPoint(i, j))
			{
				points[k].x = i;
				points[k].y = j;
				k++;
			}
		}
	}
}

//~ Selecione aleatoriamente um ponto gerador (G)
void generateG()
{
	//~ Mudar a semente   
	srand(time(NULL));
	G = points[rand() % P];
}

//~ Gerar as chaves públicas e privadas
void generateKeys()
{
	for (int i = 0; i < N; i++)
	{
		bool valid = true;
		do
		{
			//~ Chave privada (k)
			srand(time(NULL));
			PEs[i].k = mod(rand());

			//~ Chaves Públicas
			struct coord pu = eccMultBinary(PEs[i].k, G);
			if (isValidPoint(pu.x, pu.y))
			{
				PEs[i].pu = pu;
				valid = true;
				printf("\nPE%d => kG=PU => %d[%d,%d]=[%d,%d]", i, PEs[i].k, G.x, G.y,
					   pu.x, pu.y);
			}
			else
			{
				printf("\nPE%d => kG=PU => %d[%d,%d]=[%d,%d] - INVÁLIDO!", i,
					   PEs[i].k, G.x, G.y, pu.x, pu.y);
				valid = false;
			}
		}
		while (!valid);
	}
}

void showPoints()
{
	printf("\nPoints:\n");
	for (int i = 0; i < P; i++)
		printf("[%d,%d]", points[i].x, points[i].y);
}

//~ It's show time
void main()
{
	printf("\n");

	calculatePoints();
	generateG();
	generateKeys();
	showPoints();

	printf("\n\n\tG:\t\t[%d,%d]\n", G.x, G.y);
	printf("\tPU[kgc]:\t[%d,%d]\n", PEs[PKG].pu.x, PEs[PKG].pu.y);
	printf("\tPU[Alice]:\t[%d,%d]\n", PEs[ALICE].pu.x, PEs[ALICE].pu.y);
	printf("\tPU[Bob]:\t[%d,%d]\n", PEs[BOB].pu.x, PEs[BOB].pu.y);

	printf("\n");
}
