#include "ruby.h"
#include <stdlib.h>

// Defining a space for information and references about the module to be
// stored internally
VALUE Contig = Qnil;
VALUE Transrate = Qnil;

// Prototype for the initialization method - Ruby calls this, not you
void Init_transrate();

// methods are prefixed by 'method_' here
//VALUE TestInit(VALUE, VALUE, VALUE, VALUE, VALUE);
VALUE method_composition(VALUE, VALUE);
VALUE method_base_count(VALUE,VALUE);
VALUE method_dibase_count(VALUE,VALUE);
VALUE method_kmer_count(VALUE,VALUE,VALUE);
VALUE method_longest_orf(VALUE, VALUE);

int * base_counts;
int * dibase_counts;

// The initialization method for this module
void Init_transrate() {
    Transrate = rb_define_module("Transrate");
    // VALUE rb_define_class_under(VALUE outer, const char *name, VALUE super)
    Contig = rb_define_class_under(Transrate, "Contig", rb_cObject);
    // rb_define_method(Contig, "initialize", TestInit, 2);
    rb_define_method(Contig, "composition", method_composition, 1);
    rb_define_method(Contig, "base_count", method_base_count, 1);
    rb_define_method(Contig, "dibase_count", method_dibase_count, 1);
    rb_define_method(Contig, "kmer_count", method_kmer_count, 2);
    rb_define_method(Contig, "longest_orf", method_longest_orf, 1);
}

VALUE method_composition(VALUE self, VALUE _seq) {
    int i,len, idx;
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
// either the start of the sequence or a stop codon and either the
// end of the sequence or a stop codon
VALUE method_longest_orf(VALUE self, VALUE _s) {
    int i,sl,longest=0;
    int len[6];
    char * c_str;

    sl = RSTRING_LEN(_s);
    c_str = StringValueCStr(_s);
    for (i=0;i<6;i++) {
        len[i]=0;
    }
    for (i=0;i<sl-2;i++) {
        if (c_str[i]=='T' &&
        ((c_str[i+1]=='A' && c_str[i+2]=='G') ||
        (c_str[i+1]=='A' && c_str[i+2]=='A') ||
        (c_str[i+1]=='G' && c_str[i+2]=='A'))) {
            if (len[i%3] > longest) {
                longest = len[i%3];
            }
            len[i%3]=0;
        } else {
            len[i%3]++;
        }
        if (c_str[i+2]=='A' &&
        ((c_str[i]=='C' && c_str[i+1]=='T') ||
        (c_str[i]=='T' && c_str[i+1]=='T') ||
        (c_str[i]=='T' && c_str[i+1]=='C'))) {
            if (len[3+i%3] > longest) {
                longest = len[3+i%3];
            }
            len[3+i%3]=0;
        } else {
            len[3+i%3]++;
        }
    }
    if (len[i%3] > longest) {
        longest = len[i%3];
    }
    if (len[3+i%3] > longest) {
        longest = len[3+i%3];
    }
    return INT2NUM(longest);
}