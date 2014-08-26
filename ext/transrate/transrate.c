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
VALUE method_load_bcf(VALUE, VALUE, VALUE);

VALUE method_get_len(VALUE, VALUE);
VALUE method_get_contig_name(VALUE, VALUE);
VALUE method_get_mean_effective_coverage(VALUE, VALUE);
VALUE method_get_coverage_variance(VALUE, VALUE);
VALUE method_get_effective_coverage_variance(VALUE, VALUE);
VALUE method_get_uncovered_bases(VALUE, VALUE);
VALUE method_get_low_mapq_bases(VALUE, VALUE);
VALUE method_get_total_mapq(VALUE, VALUE);
VALUE method_get_total_coverage(VALUE self, VALUE _i);
VALUE method_free_contigs(VALUE self);

void calculate_metrics(int, int, int *, int *);

int * base_counts;
int * dibase_counts;

typedef struct {
    double mean_effective_coverage;
    double coverage_variance;
    double effective_coverage_variance;
    double mapq_geo_mean;
    int len;
    int total_coverage;
    int n_uncovered_bases;
    int n_low_mapq_bases;
    int total_mapq;
    char * cname;
} CONTIGINFO;

CONTIGINFO * contigs;

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
  rb_define_method(ReadMetrics, "load_bcf", method_load_bcf, 2);
  rb_define_method(ReadMetrics, "get_len",
                   method_get_len, 1);
  rb_define_method(ReadMetrics, "get_contig_name",
                   method_get_contig_name, 1);
  rb_define_method(ReadMetrics, "get_mean_effective_coverage",
                   method_get_mean_effective_coverage, 1);
  rb_define_method(ReadMetrics, "get_coverage_variance",
                   method_get_coverage_variance, 1);
  rb_define_method(ReadMetrics, "get_effective_coverage_variance",
                   method_get_effective_coverage_variance, 1);
  rb_define_method(ReadMetrics, "get_uncovered_bases",
                   method_get_uncovered_bases, 1);
  rb_define_method(ReadMetrics, "get_low_mapq_bases",
                   method_get_low_mapq_bases, 1);
  rb_define_method(ReadMetrics, "get_total_mapq",
                   method_get_total_mapq, 1);
  rb_define_method(ReadMetrics, "get_total_coverage",
                   method_get_total_coverage, 1);
  rb_define_method(ReadMetrics, "free_contigs",
                   method_free_contigs, 0);
}

VALUE method_load_bcf(VALUE self, VALUE _filename, VALUE _size) {
  FILE *fh;
  int size,pos;
  char * filename;
  size_t len;
  char * line = NULL;
  char * contig_name = NULL;
  char * previous_name = NULL;
  int * coverage_array;
  int * mapq_array;
  int i,c,d,start,end,col,p,cov,mapq,read_length;
  int num;
  int mem, prev;
  filename = StringValueCStr(_filename);
  size = NUM2INT(_size);
  contigs = malloc(size * sizeof(CONTIGINFO));
  for (i=0;i<size;i++) {
    CONTIGINFO s;
    s.len = -1;
    s.mean_effective_coverage = 0;
    s.coverage_variance = 0;
    s.effective_coverage_variance = 0;
    s.mapq_geo_mean = 0;
    s.total_coverage = 0;
    s.n_uncovered_bases = -1;
    s.n_low_mapq_bases = 0;
    s.total_mapq = 0;
    s.cname = "na";
    contigs[i] = s;
  }
  num = -1;
  pos = 0;
  mem = 0;
  prev = 0;
  read_length = 100;
  previous_name = "na";
  fh = fopen(filename, "r"); // open filename for reading
  cov=-1;
  mapq=-1;
  while ( getline(&line, &len, fh) != -1) {
    // header info on contigs
    if (line[0]=='#') {
      if (line[2]=='c') {
        num++;
        c = 0;
        start=0;
        end=0;
        while (line[c] != '\0') { // loop over header line
          if (line[c]=='<') {
            start = c+4;
          }
          if (line[c]==',') {
            end = c-1;
            contigs[num].cname = malloc((end-start+2) * sizeof(char));
            i=0;
            for (d=start;d<=end;d++) {
              contigs[num].cname[i] = line[d];
              i++;
            }
            contigs[num].cname[i] = '\0';
            start = c+8;
          }
          if (line[c]=='>') {
            end = c-1;
            char n[start-end+2];
            i=0;
            for (d=start; d<=end; d++) {
              n[i]=line[d];
              i++;
            }
            n[i]='\0';
            contigs[num].len = atoi(n);
            contigs[num].n_uncovered_bases = contigs[num].len;
          } // endif line[c]=='>'
          c++;
        } // end of while loop
      }
      if (line[1]=='C') { // last line of header. set counter back to start
        num=0;
      }
    } else { // line doesn't start with a #

      // column info
      // when found new contig
      //   then do calculations over array
      //   make new empty array of len l (get length from struct)
      //   iterate over lines and fill array with coverge and mapq data

      c = 0;
      col = 0;
      start = 0;
      char * c_string;
      while (line[c] != '\0') { // loop through line_
        if (line[c] == '\t') {
          end = c;
          if (col==0) { // name
            contig_name = malloc((end-start+1) * sizeof(char));
            i=0;
            for(d = start; d < end; d++) {
              contig_name[i]=line[d];
              i++;
            }
            contig_name[i]='\0';
            mem=i+1;
          }
          if (col==1) { // position
            c_string = malloc((end-start+1) * sizeof(char));
            i=0;
            for(d = start; d < end; d++) {
              c_string[i]=line[d];
              i++;
            }
            c_string[i]='\0';
            pos = atoi(c_string) - 1;
            free(c_string);
          }
          start = c+1;
          col++;
        }

        if (col==7) { // info, semicolon separated string
          if (line[c]==';') {
            end = c;
            if (line[start]=='D' && line[start+1]=='P' && line[start+2]=='=') {
              c_string = malloc((end-start-2) * sizeof(char));
              i=0;
              for(d = start+3; d < end; d++) {
                c_string[i]=line[d];
                i++;
              }
              c_string[i]='\0';
              cov = atoi(c_string);
              free(c_string);
            }
            if (line[start]=='M' && line[start+1]=='Q' && line[start+2]=='=') {
              c_string = malloc((end-start-2) * sizeof(char));
              i=0;
              for(d = start+3; d < end; d++) {
                c_string[i]=line[d];
                i++;
              }
              c_string[i]='\0';
              mapq = atoi(c_string);
              free(c_string);
            }
            start = c + 1;
          } // endif line[c]==';'
        } // endif col==7

        c++;
      } // end of while loop through line
      if (strcmp(previous_name, contig_name)!=0) {
        if (strcmp(previous_name, "na")==0) {
          // first line of file
        } else {
          calculate_metrics(num, read_length, coverage_array, mapq_array);
          free(coverage_array);
          free(mapq_array);
        }
        // scan through contigs array to find struct with cname==contig_name
        while (strcmp(contig_name, contigs[num].cname)!=0) {
          num++;
        }
        if (prev>0) {
          free(previous_name);
        }
        previous_name = malloc(mem * sizeof(char));
        strcpy(previous_name,contig_name);
        prev = 1;

        // get length of contig from struct array
        p = contigs[num].len;
        // create new empty array for coverage and mapq
        coverage_array = malloc(p * sizeof(int));
        mapq_array = malloc(p * sizeof(int));
        for(i=0;i<p;i++) {
          coverage_array[i]=0;
          mapq_array[i]=-1;
        }
      }
      // add cov and mapq from this line to the array
      if (cov >= 0) {
        coverage_array[pos] = cov;
      }
      if (mapq >= 0) {
        mapq_array[pos] = mapq;
      }
      if (mem > 0) {
        free(contig_name);
        mem = 0;
      }

    } // endif
  } // end of while loop through file
  free(line);
  fclose(fh);
  calculate_metrics(num, read_length, coverage_array, mapq_array);
  free(coverage_array);
  free(mapq_array);
  return INT2NUM(0);
}

void calculate_metrics(int num, int read_length, int * coverage_array,
                       int * mapq_array) {

  int c, e, i, l, m, n_uncovered_bases,n_low_mapq_bases;
  c = 0; e = 0; m = 0; n_uncovered_bases = 0;
  l = contigs[num].len;
  double mean, eff_mean, var, eff_var, mapq_geo_mean;
  var = 0; eff_var = 0;
  //coverage
  for(i=0;i<l;i++) {
    c += coverage_array[i];
    if (i >= read_length && i < contigs[num].len-read_length) {
      e += coverage_array[i];
      if (coverage_array[i]==0) {
        n_uncovered_bases++;
      }
    }
  }
  mean = c / (double)l;
  if (l>2*read_length) {
    eff_mean =  e / (double)(l-2*read_length);
  } else {
    eff_mean = 0;
  }
  // variance
  for(i=0;i<l;i++) {
    var += pow(coverage_array[i] - mean,2);
    if (i >= read_length && i < contigs[num].len-read_length) {
      eff_var += pow(coverage_array[i] - eff_mean,2);
    }
  }

  var = var / (double)l;
  if (l>2*read_length) {
    eff_var = eff_var / (double)(l-2*read_length);
  } else {
    eff_var = 0;
  }
  //mapq
  n_low_mapq_bases = 0;
  mapq_geo_mean=0;
  for(i=0;i<l;i++) {
    m += mapq_array[i];
    if (i >= read_length && i < contigs[num].len-read_length) {
      mapq_geo_mean += log(mapq_array[i]);
      if (mapq_array[i]>=0 && mapq_array[i]<5) {
        n_low_mapq_bases++;
      }
    }
  }
  //geometric mean
  mapq_geo_mean /= l;
  mapq_geo_mean = exp(mapq_geo_mean);

  // set values in contiginfo array
  contigs[num].n_uncovered_bases = n_uncovered_bases;
  contigs[num].n_low_mapq_bases = n_low_mapq_bases;
  contigs[num].total_coverage = c;
  contigs[num].mean_effective_coverage = eff_mean;
  contigs[num].coverage_variance = var;
  contigs[num].effective_coverage_variance = eff_var;
  contigs[num].total_mapq = m;
  contigs[num].mapq_geo_mean = mapq_geo_mean;

}

VALUE method_get_len(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return INT2NUM(contigs[i].len);
}

VALUE method_get_contig_name(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return rb_str_new2(contigs[i].cname);
}

VALUE method_get_mean_effective_coverage(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return rb_float_new(contigs[i].mean_effective_coverage);
}

VALUE method_get_total_coverage(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return INT2NUM(contigs[i].total_coverage);
}

VALUE method_get_coverage_variance(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return rb_float_new(contigs[i].coverage_variance);
}

VALUE method_get_effective_coverage_variance(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return rb_float_new(contigs[i].effective_coverage_variance);
}

VALUE method_get_uncovered_bases(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return INT2NUM(contigs[i].n_uncovered_bases);
}

VALUE method_get_low_mapq_bases(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return INT2NUM(contigs[i].n_low_mapq_bases);
}

VALUE method_get_total_mapq(VALUE self, VALUE _i) {
  int i;
  i = NUM2INT(_i);
  return rb_float_new(contigs[i].total_mapq);
}

VALUE method_free_contigs(VALUE self) {
  free(contigs);
  return INT2NUM(0);
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

