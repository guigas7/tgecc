struct coord
{
	int x;
	int y;
};

// Tamanho do P. Também ditará a ordem da curva.
#define P 23

// Parâmetro A da curva.
#define A 1

// Parâmetro B da curva.
#define B 0

// Número de PEs no sistema local
#define N 3

//~ Operação modular
int mod(int a)
{
	int b = (a < 0) ? (a += P) % P : a % P;	
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
	
	// Slope
	int s = mod(mod(p.y - q.y) * mod(multInv(mod(p.x - q.y))));
	// X
	r.x = mod(mod(mod(s * s) - p.x) - q.x);
	// Y
	r.y = mod(mod(s * mod(p.x - r.x)) - p.y);

	//~ printf("\n[%.2d,%.2d]+[%.2d,%.2d] -> [%.2d,%.2d]", p.x, p.y, q.x, q.y, r.x, r.y);

	return r;
}
//~ ECCDBL r=2p
struct coord eccDbl(struct coord p)
{
	struct coord r;
	// Propriedades do grupo multiplicativo. 

	if (p.y == -p.y || p.y == 0)
	{
		r.x = r.y = 0;
	}
	else if(p.y != 0)
	{
		// Slope
		int s = (mod(mod(3 * mod(p.x * p.x)) + A)) * multInv(mod(2 * p.y));
		// X
		r.x = mod(mod(s * s) - mod(2 * p.x));
		// Y
		r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
		//~ r.y = mod(p.y + mod(s * mod(r.x - p.x)));
	}
	//~ printf("\n[%.2d,%.2d]+[%.2d,%.2d] -> [%.2d,%.2d]", p.x, p.y, p.x, p.y, r.x, r.y);
	return r;
}
