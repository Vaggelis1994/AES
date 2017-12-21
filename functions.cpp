#include <lela/solutions/echelon-form.h>

using namespace LELA;

template <class Ring, class Matrix>
void ComputeRowEchelonForm (const Ring &R, Matrix &A)
{
	Context<Ring> ctx (R);
	EchelonForm<Ring> EF (ctx);

	echelonize (A); // A now replaced by its row echelon form
}