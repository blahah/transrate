#include "ruby.h"
#include <stdlib.h>
#include <math.h>

// Defining a space for information and references about the module to be
// stored internally
VALUE Contig = Qnil;
VALUE ReadMetrics = Qnil;
VALUE Transrate = Qnil;

// Prototype for the initialization method - Ruby calls this, not you
void Init_transrate();

// methods are prefixed by 'method_' here
// contig
VALUE method_composition(VALUE, VALUE);
VALUE method_base_count(VALUE,VALUE);
VALUE method_dibase_count(VALUE,VALUE);
VALUE method_kmer_count(VALUE,VALUE,VALUE);
VALUE method_longest_orf(VALUE, VALUE);
// read_metrics

int * base_counts;
int * dibase_counts;

// The initialization method for this module
void Init_transrate() {
  Transrate = rb_define_module("Transrate");
  Contig = rb_define_class_under(Transrate, "Contig", rb_cObject);
  ReadMetrics = rb_define_class_under(Transrate, "ReadMetrics", rb_cObject);
  // contig
  rb_define_method(Contig, "composition", method_composition, 1);
  rb_define_method(Contig, "base_count", method_base_count, 1);
  rb_define_method(Contig, "dibase_count", method_dibase_count, 1);
  rb_define_method(Contig, "kmer_count", method_kmer_count, 2);
  rb_define_method(Contig, "longest_orf", method_longest_orf, 1);
  // ReadMetrics
}

VALUE method_composition(VALUE self, VALUE _seq) {
  int i, len, idx;
  char * seq;
  char base;
  char prevbase;
  seq = StringValueCStr(_seq);
  len = RSTRING_LEN(_seq);
  base_counts = malloc(5 * sizeof(int));
  dibase_counts = malloc(25 * sizeof(int));

  for (i=0; i < 5; i++) {
    base_counts[i]=0;
  }
  for (i=0; i < 25; i++) {
    dibase_counts[i]=0;
  }
  for (i=0; i < len; i++) {
    base = seq[i];
    if (base > 90) {
      base -= 32;
    }
    switch (base) {
      case 'A': {
        idx=0;
        break;
      }
      case 'C': {
        idx=1;
        break;
      }
      case 'G': {
        idx=2;
        break;
      }
      case 'T': {
        idx=3;
        break;
      }
      default: {
        idx=4;
        break;
      }
    }
    base_counts[idx]++;

    if (i > 0) {
      prevbase = seq[i-1];
      if (prevbase > 90) {
        prevbase -= 32;
      }
      switch (prevbase) {
        case 'A': {
          idx=idx;
          break;
        }
        case 'C': {
          idx=idx+5;
          break;
        }
        case 'G': {
          idx=idx+10;
          break;
        }
        case 'T': {
          idx=idx+15;
          break;
        }
        default: {
          idx=idx+20;
          break;
        }
      }
      dibase_counts[idx]++;
    }
  }
  return INT2NUM(0);
}

VALUE method_dibase_count(VALUE self, VALUE idx) {
  return INT2NUM(dibase_counts[NUM2INT(idx)]);
}

VALUE method_base_count(VALUE self, VALUE idx) {
  return INT2NUM(base_counts[NUM2INT(idx)]);
}

VALUE method_kmer_count(VALUE self, VALUE _k, VALUE _s) {
  int n, i, start, k, len, h, size = 0;
  char * c_str;
  char base;
  len = RSTRING_LEN(_s);
  c_str = StringValueCStr(_s);
  k = NUM2INT(_k);
  size = 1;
  for(h=0;h<k;h++) {
      size *= 4;
  }
  short set[size];
  for(start=0;start<size;start++) {
      set[start]=0;
  }
  for(start=0; start<len-k+1; start++) {
    i = 0;
    h = 0;
    n = 0;
    for(i = start; i < start+k; i++) {
      base = c_str[i];
      if (base > 90) {
        base -= 32;
      }
      switch (base) {
        case 'A': {
          h = h << 2;
          h += 0;
          break;
        }
        case 'C': {
          h = h << 2;
          h += 1;
          break;
        }
        case 'G': {
          h = h << 2;
          h += 2;
          break;
        }
        case 'T': {
          h = h << 2;
          h += 3;
          break;
        }
        default: {
          n++;
          break;
        }
      }
    }
    if (n==0) {
      set[h] += 1;
    }
  }
  i = 0; // count how many in array are set //
  for(start = 0; start < size; start++) {
    if (set[start]>0) {
      i++;
    }
  }
  return INT2NUM(i);
}

// takes in a string and calculates the longest open reading frame
// in any of the 6 frames
// an open reading frame is defined as the number of bases between
// either the start of the sequence or a start codon and either the
// end of the sequence or a stop codon

VALUE method_longest_orf(VALUE self, VALUE _str) {
  int i,sl,longest=0;
  int len[3];
  char * str;
  sl = RSTRING_LEN(_str);
  str = StringValueCStr(_str);
  for (i=0;i<3;i++) {
    len[i]=0;
  }
  for(i=0;i<sl-2;i++) {
    if (str[i]=='A' && str[i+1]=='T' && str[i+2]=='G') { //Methionine
      if (len[i%3]>=0) {
        len[i%3]++;
      } else {
        len[i%3]=1;
      }
    } else {
      if (str[i]=='T' &&
        ((str[i+1]=='A' && str[i+2]=='G') ||   //amber
        (str[i+1]=='A' && str[i+2]=='A') ||    //ochre   stops
        (str[i+1]=='G' && str[i+2]=='A'))) {   //umber
        if (len[i%3]>longest) {
          longest = len[i%3];
        }
        len[i%3]=-1;
      } else { // any other codon
        if (len[i%3]>=0) {
          len[i%3]++;
        }
      }
    }
  }
  for(i=0;i<3;i++) {
    if (len[i%3] > longest) {
      longest = len[i%3];
    }
  }
  for (i=0;i<3;i++) {
    len[i]=0;
  }
  for(i=sl-1;i>=2;i--) {
    if (str[i]=='T' && str[i-1]=='A' && str[i-2]=='C') { //Methionine
      if (len[i%3]>=0) {
        len[i%3]++;
      } else {
        len[i%3]=1;
      }
    } else {
      if (str[i]=='A' &&
        ((str[i-1]=='T' && str[i-2]=='C') ||   //amber
        (str[i-1]=='T' && str[i-2]=='T') ||    //ochre   stops
        (str[i-1]=='C' && str[i-2]=='T'))) {   //umber
        if (len[i%3]>longest) {
          longest = len[i%3];
        }
        len[i%3]=-1;
      } else { // any other codon
        if (len[i%3]>=0) {
          len[i%3]++;
        }
      }
    }
  }
  for(i=0;i<3;i++) {
    if (len[i%3] > longest) {
      longest = len[i%3];
    }
  }
  return INT2NUM(longest);
}

