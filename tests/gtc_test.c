/* gtc.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) genetic toolchain calculator (GTC).

This Test class defines cases to verify that we don't break the excepted behaviours in the future upon changes.

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#include "../gtc.h" /* Genetic Toolchain Calculator */
#include "test.h"   /* Simple Testing framework */

int gtc_test_strcmp(char *a, char *b)
{
  while (*a && (*a == *b))
  {
    a++;
    b++;
  }
  return *(unsigned char *)a - *(unsigned char *)b;
}

void gtc_test_dna_to_rna(void)
{
  unsigned char dna[] = "ATGCGTATTA";
  unsigned char rna[32];
  gtc_dna_to_rna(dna, rna, sizeof(rna));
  assert(gtc_test_strcmp((char *)rna, "AUGCGUAUUA") == 0);
}

void gtc_test_base_index_rna(void)
{
  assert(gtc_base_index_rna('A') == 0);
  assert(gtc_base_index_rna('C') == 1);
  assert(gtc_base_index_rna('G') == 2);
  assert(gtc_base_index_rna('U') == 3);
  assert(gtc_base_index_rna('T') == 0); /* not RNA */
}

void gtc_test_codon_index_rna(void)
{
  /* Codon AUG = A (0), U (3), G (2)      */
  /* Index = (0 << 4) | (3 << 2) | 2 = 14 */
  assert(gtc_codon_index_rna('A', 'U', 'G') == 14);
  /* Codon CCC = 1, 1, 1 → (1<<4) | (1<<2) | 1 = 21 */
  assert(gtc_codon_index_rna('C', 'C', 'C') == 21);
}

void gtc_test_rna_to_protein(void)
{
  /* RNA: AUG GCC UUU UAA (M, A, F, stop) */
  unsigned char rna[] = "AUGGCCUUUUAA";
  unsigned char protein[32];
  gtc_rna_to_protein(rna, protein, sizeof(protein));
  assert(gtc_test_strcmp((char *)protein, "MAF") == 0);
}

void gtc_test_dna_is_valid(void)
{
  unsigned char valid1[] = "ACGTACGT";
  unsigned char valid2[] = "T";
  unsigned char empty[] = "";
  unsigned char invalid1[] = "AXGT"; /* X not valid      */
  unsigned char invalid2[] = "123";  /* digits not valid */

  assert(gtc_dna_is_valid(valid1) == 1);
  assert(gtc_dna_is_valid(valid2) == 1);
  assert(gtc_dna_is_valid(empty) == 1);
  assert(gtc_dna_is_valid(invalid1) == 0);
  assert(gtc_dna_is_valid(invalid2) == 0);
}

void gtc_test_protein_is_valid(void)
{
  unsigned char valid1[] = "ACDEFGHIKLMNPQRSTVWY";
  unsigned char valid2[] = "MAF";
  unsigned char empty[] = "";
  unsigned char invalid1[] = "MAFX"; /* X not valid    */
  unsigned char invalid2[] = "123";  /* digits invalid */

  assert(gtc_protein_is_valid(valid1) == 1);
  assert(gtc_protein_is_valid(valid2) == 1);
  assert(gtc_protein_is_valid(empty) == 1);
  assert(gtc_protein_is_valid(invalid1) == 0);
  assert(gtc_protein_is_valid(invalid2) == 0);
}

void gtc_test_protein_weight(void)
{
  unsigned char protein1[] = "MAF"; /* M = 149, A = 89, F = 165 → total = 403 */
  unsigned char protein2[] = "GTC"; /* G = 75, T = 119, C = 121 → total = 315 */
  unsigned char protein3[] = "";    /* empty string = 0                       */
  unsigned char invalid[] = "XQZ";  /* X and Z invalid = weight = 0           */

  assert(gtc_protein_weight(protein1) == (149 + 89 + 165)); /* 403 */
  assert(gtc_protein_weight(protein2) == (75 + 119 + 121)); /* 315 */
  assert(gtc_protein_weight(protein3) == 0);
  assert(gtc_protein_weight(invalid) == 0);
}

void gtc_test_reverse_complement(void)
{
  unsigned char dna1[] = "ATCG";
  unsigned char out1[5];
  unsigned char dna2[] = "GATTACA";
  unsigned char out2[8];
  unsigned char dna3[] = "";
  unsigned char out3[1];

  gtc_reverse_complement(dna1, out1, 5);
  assert(gtc_test_strcmp((char *)out1, "CGAT") == 0);

  gtc_reverse_complement(dna2, out2, 8);
  assert(gtc_test_strcmp((char *)out2, "TGTAATC") == 0);

  gtc_reverse_complement(dna3, out3, 1);
  assert(gtc_test_strcmp((char *)out3, "") == 0);
}

void gtc_test_gc_content_percent(void)
{
  unsigned char *seq1 = (unsigned char *)"";
  unsigned char *seq2 = (unsigned char *)"ATATAT";
  unsigned char *seq3 = (unsigned char *)"GCGCGC";
  unsigned char *seq4 = (unsigned char *)"AGCT";
  unsigned char *seq5 = (unsigned char *)"GGATCC";

  assert(gtc_gc_content_percent(seq1) == 0);   /* Empty */
  assert(gtc_gc_content_percent(seq2) == 0);   /* 0% GC */
  assert(gtc_gc_content_percent(seq3) == 100); /* 100% GC */
  assert(gtc_gc_content_percent(seq4) == 50);  /* 2 of 4 GC */
  assert(gtc_gc_content_percent(seq5) == 66);  /* 4 of 6 GC */
}

void gtc_test_count_codon_occurrences(void)
{
  unsigned char rna1[] = "AUGGCUAUGGCU";
  unsigned char rna2[] = "AUGAUGAUG";
  unsigned char rna3[] = "";
  unsigned char codon1[3] = {'A', 'U', 'G'};
  unsigned char codon2[3] = {'G', 'C', 'U'};

  assert(gtc_count_codon_occurrences(rna1, codon1) == 2);
  assert(gtc_count_codon_occurrences(rna1, codon2) == 2);
  assert(gtc_count_codon_occurrences(rna2, codon1) == 3);
  assert(gtc_count_codon_occurrences(rna3, codon1) == 0);
}

int main(void)
{
  gtc_test_dna_to_rna();
  gtc_test_base_index_rna();
  gtc_test_codon_index_rna();
  gtc_test_rna_to_protein();
  gtc_test_dna_is_valid();
  gtc_test_protein_is_valid();
  gtc_test_protein_weight();
  gtc_test_reverse_complement();
  gtc_test_gc_content_percent();
  gtc_test_count_codon_occurrences();

  return 0;
}

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
