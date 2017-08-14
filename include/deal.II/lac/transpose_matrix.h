// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii__transpose_matrix_h
#define dealii__transpose_matrix_h


#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/pointer_matrix.h>

DEAL_II_NAMESPACE_OPEN


/**
 * The transpose of a given matrix.  This auxiliary class swaps the effect of
 * vmult() and Tvmult() as well as vmult_add() and Tvmult_add().
 *
 * The implementation is analogous to the class PointerMatrix.
 *
 * @note The transposed matrix is never actually assembled. Instead, only the
 * matrix vector multiplication is performed in a transposed way.
 *
 * @deprecated If deal.II was configured with C++11 support, use the
 * LinearOperator class instead, see the module on
 * @ref LAOperators "linear operators"
 * for further details.
 *
 * @ingroup Matrix2
 * @author Guido Kanschat, 2006
 */
template<typename MatrixType, typename VectorType>
class
  TransposeMatrix : public PointerMatrixBase<VectorType>
{
public:
  /**
   * Constructor.  The pointer in the argument is stored in this class. As
   * usual, the lifetime of <tt>*M</tt> must be longer than the one of the
   * PointerMatrix.
   *
   * If <tt>M</tt> is zero, no matrix is stored.
   */
  TransposeMatrix (const MatrixType *M=0);

  /**
   * Constructor. The name argument is used to identify the SmartPointer for
   * this object.
   */
  TransposeMatrix(const char *name);

  /**
   * Constructor. <tt>M</tt> points to a matrix which must live longer than
   * the TransposeMatrix. The name argument is used to identify the
   * SmartPointer for this object.
   */
  TransposeMatrix(const MatrixType *M,
                  const char       *name);

  // Use doc from base class
  virtual void clear();

  /**
   * Return whether the object is empty.
   */
  bool empty () const;

  /**
   * Assign a new matrix pointer. Deletes the old pointer and releases its
   * matrix.
   * @see SmartPointer
   */
  const TransposeMatrix &operator= (const MatrixType *M);

  /**
   * Matrix-vector product.
   */
  virtual void vmult (VectorType       &dst,
                      const VectorType &src) const;

  /**
   * Transposed matrix-vector product.
   */
  virtual void Tvmult (VectorType       &dst,
                       const VectorType &src) const;

  /**
   * Matrix-vector product, adding to <tt>dst</tt>.
   */
  virtual void vmult_add (VectorType       &dst,
                          const VectorType &src) const;

  /**
   * Transposed matrix-vector product, adding to <tt>dst</tt>.
   */
  virtual void Tvmult_add (VectorType       &dst,
                           const VectorType &src) const;

private:
  /**
   * The pointer to the actual matrix.
   */
  SmartPointer<const MatrixType,TransposeMatrix<MatrixType,VectorType> > m;
};


//----------------------------------------------------------------------//


template<typename MatrixType, typename VectorType>
TransposeMatrix<MatrixType, VectorType>::TransposeMatrix (const MatrixType *M)
  : m(M)
{}


template<typename MatrixType, typename VectorType>
TransposeMatrix<MatrixType, VectorType>::TransposeMatrix (const char *name)
  : m(0, name)
{}


template<typename MatrixType, typename VectorType>
TransposeMatrix<MatrixType, VectorType>::TransposeMatrix (const MatrixType *M,
                                                          const char       *name)
  : m(M, name)
{}


template<typename MatrixType, typename VectorType>
inline void
TransposeMatrix<MatrixType, VectorType>::clear ()
{
  m = 0;
}


template<typename MatrixType, typename VectorType>
inline const TransposeMatrix<MatrixType, VectorType> &
TransposeMatrix<MatrixType, VectorType>::operator= (const MatrixType *M)
{
  m = M;
  return *this;
}


template<typename MatrixType, typename VectorType>
inline bool
TransposeMatrix<MatrixType, VectorType>::empty () const
{
  if (m == 0)
    return true;
  return m->empty();
}

template<typename MatrixType, typename VectorType>
inline void
TransposeMatrix<MatrixType, VectorType>::vmult (VectorType       &dst,
                                                const VectorType &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult (dst, src);
}


template<typename MatrixType, typename VectorType>
inline void
TransposeMatrix<MatrixType, VectorType>::Tvmult (VectorType       &dst,
                                                 const VectorType &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult (dst, src);
}


template<typename MatrixType, typename VectorType>
inline void
TransposeMatrix<MatrixType, VectorType>::vmult_add (VectorType       &dst,
                                                    const VectorType &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->Tvmult_add (dst, src);
}


template<typename MatrixType, typename VectorType>
inline void
TransposeMatrix<MatrixType, VectorType>::Tvmult_add (VectorType       &dst,
                                                     const VectorType &src) const
{
  Assert (m != 0, ExcNotInitialized());
  m->vmult_add (dst, src);
}


DEAL_II_NAMESPACE_CLOSE

#endif