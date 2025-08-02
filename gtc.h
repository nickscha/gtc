/* gtc.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) genetic toolchain calculator (GTC).

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#ifndef GTC_H
#define GTC_H

/* #############################################################################
 * # COMPILER SETTINGS
 * #############################################################################
 */
/* Check if using C99 or later (inline is supported) */
#if __STDC_VERSION__ >= 199901L
#define GTC_INLINE inline
#define GTC_API extern
#elif defined(__GNUC__) || defined(__clang__)
#define GTC_INLINE __inline__
#define GTC_API static
#elif defined(_MSC_VER)
#define GTC_INLINE __inline
#define GTC_API static
#else
#define GTC_INLINE
#define GTC_API static
#endif

/* === Check if DNA sequence is valid (only A, C, G, T allowed) === */
GTC_API GTC_INLINE int gtc_dna_is_valid(unsigned char *dna)
{
  unsigned short i = 0;
  unsigned char b;

  while ((b = dna[i]) != 0)
  {
    if (b != 'A' && b != 'C' && b != 'G' && b != 'T')
    {
      return 0; /* invalid character found */
    }
    i++;
  }

  return 1; /* valid */
}

/* === DNA -> RNA Transcription === */
GTC_API GTC_INLINE void gtc_dna_to_rna(unsigned char *dna, unsigned char *rna, unsigned short max_len)
{
  unsigned short i;

  for (i = 0; i < max_len - 1; ++i)
  {
    unsigned char b = dna[i];

    if (b == 0)
    {
      rna[i] = 0;
      return;
    }

    rna[i] = (b == 'T') ? 'U' : b;
  }

  rna[max_len - 1] = 0;
}

/* === RNA -> Protein Translation === */
static unsigned char gtc_codon_table[64] = {
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'};

GTC_API GTC_INLINE unsigned char gtc_base_index_rna(unsigned char b)
{
  switch (b)
  {
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'U':
    return 3;
  default:
    break;
  }

  return 0;
}

GTC_API GTC_INLINE unsigned char gtc_codon_index_rna(unsigned char b1, unsigned char b2, unsigned char b3)
{
  return (unsigned char)((gtc_base_index_rna(b1) << 4) | (gtc_base_index_rna(b2) << 2) | gtc_base_index_rna(b3));
}

GTC_API GTC_INLINE void gtc_rna_to_protein(unsigned char *rna, unsigned char *protein, unsigned short max_len)
{
  unsigned short i = 0, j = 0;

  while (rna[i] && rna[i + 1] && rna[i + 2] && j < max_len - 1)
  {
    unsigned char idx = gtc_codon_index_rna(rna[i], rna[i + 1], rna[i + 2]);
    unsigned char amino_acid = gtc_codon_table[idx];

    if (amino_acid == '*')
    {
      break;
    }

    protein[j++] = amino_acid;
    i += 3;
  }

  protein[j] = 0;
}

/* map from amino acid to weight (scaled by 10) */
static unsigned short gtc_protein_weights[26] = {
    /* A B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z */
    89, 0, 121, 133, 147, 165, 75, 155, 131, 0, 146, 131, 149, 132, 0, 115, 146, 174, 105, 119, 0, 117, 204, 0, 181, 0};

/* === Molecular weights lookup for 20 standard amino acids (units: Daltons / 10) === */
GTC_API GTC_INLINE unsigned short gtc_protein_weight(unsigned char *protein)
{
  unsigned short sum = 0;
  unsigned short i = 0;
  unsigned char amino_acid;

  while ((amino_acid = protein[i]) != 0)
  {
    unsigned short weight;

    if (amino_acid < 'A' || amino_acid > 'Z')
    {
      return 0; /* invalid character */
    }

    weight = gtc_protein_weights[amino_acid - 'A'];

    if (weight == 0)
    {
      return 0; /* invalid amino acid */
    }

    sum += weight;
    i++;
  }

  return sum;
}

/* === Check if protein sequence is valid (only 20 standard amino acids allowed) === */
GTC_API GTC_INLINE int gtc_protein_is_valid(unsigned char *protein)
{
  unsigned short i = 0;
  unsigned char amino_acid;

  /* Valid amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y */
  while ((amino_acid = protein[i]) != 0)
  {
    switch (amino_acid)
    {
    case 'A':
    case 'C':
    case 'D':
    case 'E':
    case 'F':
    case 'G':
    case 'H':
    case 'I':
    case 'K':
    case 'L':
    case 'M':
    case 'N':
    case 'P':
    case 'Q':
    case 'R':
    case 'S':
    case 'T':
    case 'V':
    case 'W':
    case 'Y':
      break; /* valid */
    default:
      return 0; /* invalid */
    }

    i++;
  }

  return 1;
}

GTC_API GTC_INLINE void gtc_reverse_complement(unsigned char *dna, unsigned char *revcomp, unsigned short max_len)
{
  unsigned short len = 0;

  while (dna[len] && len < max_len - 1)
  {
    len++;
  }

  if (len >= max_len)
  {
    revcomp[0] = 0;
    return;
  }

  while (len > 0)
  {
    len--;
    switch (dna[len])
    {
    case 'A':
      revcomp[max_len - len - 2] = 'T';
      break;
    case 'T':
      revcomp[max_len - len - 2] = 'A';
      break;
    case 'C':
      revcomp[max_len - len - 2] = 'G';
      break;
    case 'G':
      revcomp[max_len - len - 2] = 'C';
      break;
    default:
      revcomp[max_len - len - 2] = 'N';
      break;
    }
  }

  revcomp[max_len - 1] = 0;
}

GTC_API GTC_INLINE unsigned short gtc_gc_content_percent(unsigned char *dna)
{
  unsigned short gc = 0, total = 0;
  unsigned char b;

  while ((b = dna[total]) != 0)
  {
    if (b == 'G' || b == 'C')
    {
      gc++;
    }
    total++;
  }

  if (total == 0)
  {
    return 0;
  }

  return (unsigned short)((gc * 100) / total);
}

GTC_API GTC_INLINE unsigned short gtc_count_codon_occurrences(unsigned char *rna, unsigned char codon[3])
{
  unsigned short count = 0, i = 0;

  while (rna[i] && rna[i + 1] && rna[i + 2])
  {
    if (rna[i] == codon[0] && rna[i + 1] == codon[1] && rna[i + 2] == codon[2])
    {
      count++;
    }
    i += 3;
  }

  return count;
}

#endif /* GTC_H */

/*
   ------------------------------------------------------------------------------
   This software is available under 2 licenses -- choose whichever you prefer.
   ------------------------------------------------------------------------------
   ALTERNATIVE A - MIT License
   Copyright (c) 2025 nickscha
   Permission is hereby granted, free of charge, to any person obtaining a copy of
   this software and associated documentation files (the "Software"), to deal in
   the Software without restriction, including without limitation the rights to
   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is furnished to do
   so, subject to the following conditions:
   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   ------------------------------------------------------------------------------
   ALTERNATIVE B - Public Domain (www.unlicense.org)
   This is free and unencumbered software released into the public domain.
   Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
   software, either in source code form or as a compiled binary, for any purpose,
   commercial or non-commercial, and by any means.
   In jurisdictions that recognize copyright laws, the author or authors of this
   software dedicate any and all copyright interest in the software to the public
   domain. We make this dedication for the benefit of the public at large and to
   the detriment of our heirs and successors. We intend this dedication to be an
   overt act of relinquishment in perpetuity of all present and future rights to
   this software under copyright law.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
   ------------------------------------------------------------------------------
*/
