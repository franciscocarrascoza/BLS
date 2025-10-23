/*
  bls_231025.c — Bravais Lattice Sampling (BLS) with seed-triggered refinement + benchmark
  Date: 2025-10-23
  Summary:
    - Pre-filter via Bravais lattice, then refine with DFS/Skip-DFS only when a seed hits.
    - Lattices: cubic, hexagonal, triclinic. Centering: P, F, I.
    - Scale so the shortest site distance equals dNN. If dNN=0, auto from radii with factor alpha.
    - Benchmark mode sweeps N=100..1000 and prints a table of time and memory.

  CLI defaults:
    --nx 256 --ny 256 --nz 256
    --prob 0.02 --seed 123456789
    --lattice cubic --centering F
    --dnn 0 --alpha 0.7 --skip 3 --radii 1.0,2.0
    --hex-c-over-a 1.633
    --triclinic-a 1.0 --triclinic-b 1.2 --triclinic-c 1.4
    --triclinic-alpha 90 --triclinic-beta 100 --triclinic-gamma 110
    --bench 0 --nmin 100 --nmax 1000 --nstep 100

  Build:
    cc -O3 -march=native -Wall -Wextra -std=c11 -o bls_231025 bls_231025.c -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#ifdef __linux__
#include <sys/resource.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------- RNG ----------
static uint64_t xorshift64star(uint64_t *x){ *x ^= *x >> 12; *x ^= *x << 25; *x ^= *x >> 27; return *x * 2685821657736338717ULL; }
static double urand(uint64_t *s){ return (xorshift64star(s) >> 11) * (1.0/9007199254740992.0); }

// ---------- Params ----------
typedef struct {
  int NX, NY, NZ;
  double prob;
  uint64_t seed;
  char lattice[32];
  char centering[4];
  double dNN;
  double alpha;
  int skip;
  double *radii; int n_radii;
  double hex_c_over_a;
  double ta, tb, tc;
  double talpha_deg, tbeta_deg, tgamma_deg;
  // benchmark
  int bench;
  int nmin, nmax, nstep;
} Config;

static void set_defaults(Config *cfg){
  cfg->NX=256; cfg->NY=256; cfg->NZ=256;
  cfg->prob=0.02; cfg->seed=123456789ULL;
  strcpy(cfg->lattice,"cubic"); strcpy(cfg->centering,"F");
  cfg->dNN=0.0; cfg->alpha=0.7; cfg->skip=3;
  cfg->radii=NULL; cfg->n_radii=0;
  cfg->hex_c_over_a=1.633;
  cfg->ta=1.0; cfg->tb=1.2; cfg->tc=1.4;
  cfg->talpha_deg=90.0; cfg->tbeta_deg=100.0; cfg->tgamma_deg=110.0;
  cfg->bench=0; cfg->nmin=100; cfg->nmax=1000; cfg->nstep=100;
}

static void usage(const char* prog){
  fprintf(stderr,
    "Usage: %s [--nx N --ny N --nz N] [--prob P] [--seed S]\n"
    "           [--lattice cubic|hexagonal|triclinic] [--centering P|F|I]\n"
    "           [--dnn D] [--alpha A] [--skip K] [--radii r1,r2,...]\n"
    "           [--hex-c-over-a R]\n"
    "           [--triclinic-a A --triclinic-b B --triclinic-c C --triclinic-alpha deg --triclinic-beta deg --triclinic-gamma deg]\n"
    "           [--bench 0|1] [--nmin 100] [--nmax 1000] [--nstep 100]\n", prog);
}

// portable strdup
static char* my_strdup(const char* s){ size_t n=strlen(s)+1; char* p=(char*)malloc(n); if(p) memcpy(p,s,n); return p; }

static void parse_list_double(const char* s, double **out, int *n){
  char *dup=my_strdup(s); if(!dup){ *out=NULL; *n=0; return; }
  int cap=8, cnt=0; double *arr=(double*)malloc(cap*sizeof(double));
  for(char *tok=strtok(dup,","); tok; tok=strtok(NULL,",")){
    if(cnt==cap){ cap*=2; arr=(double*)realloc(arr,cap*sizeof(double)); }
    arr[cnt++]=atof(tok);
  }
  free(dup); *out=arr; *n=cnt;
}

static void parse_cli(int argc, char**argv, Config *cfg){
  set_defaults(cfg);
  for(int i=1;i<argc;i++){
    if(!strcmp(argv[i],"--nx") && i+1<argc) cfg->NX=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--ny") && i+1<argc) cfg->NY=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--nz") && i+1<argc) cfg->NZ=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--prob") && i+1<argc) cfg->prob=atof(argv[++i]);
    else if(!strcmp(argv[i],"--seed") && i+1<argc) cfg->seed=strtoull(argv[++i],NULL,10);
    else if(!strcmp(argv[i],"--lattice") && i+1<argc) strncpy(cfg->lattice,argv[++i],31);
    else if(!strcmp(argv[i],"--centering") && i+1<argc) strncpy(cfg->centering,argv[++i],3);
    else if(!strcmp(argv[i],"--dnn") && i+1<argc) cfg->dNN=atof(argv[++i]);
    else if(!strcmp(argv[i],"--alpha") && i+1<argc) cfg->alpha=atof(argv[++i]);
    else if(!strcmp(argv[i],"--skip") && i+1<argc) cfg->skip=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--radii") && i+1<argc) parse_list_double(argv[++i], &cfg->radii, &cfg->n_radii);
    else if(!strcmp(argv[i],"--hex-c-over-a") && i+1<argc) cfg->hex_c_over_a=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-a") && i+1<argc) cfg->ta=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-b") && i+1<argc) cfg->tb=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-c") && i+1<argc) cfg->tc=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-alpha") && i+1<argc) cfg->talpha_deg=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-beta") && i+1<argc) cfg->tbeta_deg=atof(argv[++i]);
    else if(!strcmp(argv[i],"--triclinic-gamma") && i+1<argc) cfg->tgamma_deg=atof(argv[++i]);
    else if(!strcmp(argv[i],"--bench") && i+1<argc) cfg->bench=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--nmin") && i+1<argc) cfg->nmin=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--nmax") && i+1<argc) cfg->nmax=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--nstep") && i+1<argc) cfg->nstep=atoi(argv[++i]);
    else { usage(argv[0]); exit(1); }
  }
  if(cfg->alpha<=0.0 || cfg->alpha>1.5) cfg->alpha=0.7;
  if(cfg->skip<1) cfg->skip=1;
  if(cfg->nmin<10) cfg->nmin=100;
  if(cfg->nmax<cfg->nmin) cfg->nmax=cfg->nmin;
  if(cfg->nstep<1) cfg->nstep=100;
}

// ---------- Grid ----------
typedef struct {
  int NX, NY, NZ;
  unsigned char *occ;
  unsigned char *vis;
} Grid;

static inline size_t idx3(const Grid *g,int x,int y,int z){ return ((size_t)z*g->NY + (size_t)y)*g->NX + (size_t)x; }

static void grid_init(Grid *g, int NX,int NY,int NZ){
  g->NX=NX; g->NY=NY; g->NZ=NZ;
  size_t N=(size_t)NX*NY*NZ;
  g->occ=(unsigned char*)calloc(N,1);
  g->vis=(unsigned char*)calloc(N,1);
}

static void grid_free(Grid *g){
  free(g->occ); free(g->vis);
}

static void grid_fill_random(Grid *g, double prob, uint64_t *rng){
  size_t N=(size_t)g->NX*g->NY*g->NZ;
  for(size_t i=0;i<N;i++) g->occ[i] = (urand(rng) < prob) ? 1 : 0;
}

// ---------- Linear algebra ----------
typedef struct { double m[3][3]; } Mat3; // columns are basis vectors
typedef struct { double v[3]; } Vec3;

static Vec3 v3(double x,double y,double z){ Vec3 v={{x,y,z}}; return v; }
static Vec3 v3add(Vec3 a, Vec3 b){ return v3(a.v[0]+b.v[0],a.v[1]+b.v[1],a.v[2]+b.v[2]); }
static Vec3 v3sub(Vec3 a, Vec3 b){ return v3(a.v[0]-b.v[0],a.v[1]-b.v[1],a.v[2]-b.v[2]); }
static double v3norm(Vec3 a){ return sqrt(a.v[0]*a.v[0]+a.v[1]*a.v[1]+a.v[2]*a.v[2]); }
static Vec3 M_col(const Mat3 *B,int j){ return v3(B->m[0][j],B->m[1][j],B->m[2][j]); }
static Vec3 M_mul(const Mat3 *B, int i,int j,int k){
  return v3(i*B->m[0][0] + j*B->m[0][1] + k*B->m[0][2],
            i*B->m[1][0] + j*B->m[1][1] + k*B->m[1][2],
            i*B->m[2][0] + j*B->m[2][1] + k*B->m[2][2]);
}

// primitive bases
static Mat3 basis_cubic_unit(){ Mat3 B={ .m={{1,0,0},{0,1,0},{0,0,1}} }; return B; }
static Mat3 basis_hex_unit(double c_over_a){
  Mat3 B;
  B.m[0][0]=1.0;   B.m[1][0]=0.0;             B.m[2][0]=0.0;
  B.m[0][1]=0.5;   B.m[1][1]=sqrt(3.0)/2.0;   B.m[2][1]=0.0;
  B.m[0][2]=0.0;   B.m[1][2]=0.0;             B.m[2][2]=c_over_a;
  return B;
}
static Mat3 basis_triclinic_unit(double a,double b,double c,double alpha_deg,double beta_deg,double gamma_deg){
  double alpha=alpha_deg*M_PI/180.0, beta=beta_deg*M_PI/180.0, gamma=gamma_deg*M_PI/180.0;
  Mat3 B;
  B.m[0][0]=a;   B.m[1][0]=0.0;                B.m[2][0]=0.0;
  B.m[0][1]=b*cos(gamma); B.m[1][1]=b*sin(gamma); B.m[2][1]=0.0;
  B.m[0][2]=c*cos(beta);
  double cy = (cos(alpha)*b*c - (B.m[0][1]*B.m[0][2]))/(B.m[1][1]==0.0 ? 1e-9 : B.m[1][1]);
  B.m[1][2]=cy;
  double cz_sq = c*c - B.m[0][2]*B.m[0][2] - cy*cy; if(cz_sq<1e-12) cz_sq=1e-12;
  B.m[2][2]=sqrt(cz_sq);
  return B;
}

// centering offsets (lattice coords)
typedef struct { Vec3 *ofs; int n; } Offsets;
static Offsets centering_offsets(const char *cent){
  Offsets O={0};
  if(!strcmp(cent,"P")){ O.n=1; O.ofs=(Vec3*)malloc(sizeof(Vec3)); O.ofs[0]=v3(0,0,0); }
  else if(!strcmp(cent,"F")){
    O.n=4; O.ofs=(Vec3*)malloc(sizeof(Vec3)*4);
    O.ofs[0]=v3(0,0,0); O.ofs[1]=v3(0,0.5,0.5); O.ofs[2]=v3(0.5,0,0.5); O.ofs[3]=v3(0.5,0.5,0);
  } else if(!strcmp(cent,"I")){
    O.n=2; O.ofs=(Vec3*)malloc(sizeof(Vec3)*2);
    O.ofs[0]=v3(0,0,0); O.ofs[1]=v3(0.5,0.5,0.5);
  } else { O.n=1; O.ofs=(Vec3*)malloc(sizeof(Vec3)); O.ofs[0]=v3(0,0,0); }
  return O;
}

static double shortest_unit_spacing(Mat3 B0, Offsets O){
  double dmin=1e9;
  for(int di=-1; di<=1; ++di)
    for(int dj=-1; dj<=1; ++dj)
      for(int dk=-1; dk<=1; ++dk)
        for(int p=0;p<O.n;++p)
          for(int q=0;q<O.n;++q){
            Vec3 base = M_mul(&B0, di,dj,dk);
            Vec3 delta = v3add(base, v3sub(O.ofs[p], O.ofs[q]));
            double n = v3norm(delta);
            if(n>1e-12 && n<dmin) dmin=n;
          }
  return dmin;
}

static Mat3 scale_basis(Mat3 B0, double scale){
  Mat3 B=B0;
  for(int c=0;c<3;++c){ B.m[0][c]*=scale; B.m[1][c]*=scale; B.m[2][c]*=scale; }
  return B;
}

// ---------- Refinement ----------
typedef struct { int x,y,z; } Node;

static inline bool in_bounds(const Grid *g,int x,int y,int z){
  return (unsigned)x < (unsigned)g->NX && (unsigned)y < (unsigned)g->NY && (unsigned)z < (unsigned)g->NZ;
}

static int refine_cluster(Grid *g, int sx,int sy,int sz, int skip){
  static const int dx6[6]={1,-1,0,0,0,0};
  static const int dy6[6]={0,0,1,-1,0,0};
  static const int dz6[6]={0,0,0,0,1,-1};
  size_t cap=1024, top=0; Node *stack=(Node*)malloc(cap*sizeof(Node));
  int count=0;
  stack[top++] = (Node){sx,sy,sz};
  g->vis[idx3(g,sx,sy,sz)]=1;
  while(top){
    Node cur = stack[--top];
    count++;
    for(int dir=0; dir<6; ++dir){
      for(int step=1; step<=skip; ++step){
        int nx = cur.x + step*dx6[dir];
        int ny = cur.y + step*dy6[dir];
        int nz = cur.z + step*dz6[dir];
        if(!in_bounds(g,nx,ny,nz)) break;
        size_t id = idx3(g,nx,ny,nz);
        if(!g->occ[id]) break; // stop direction on first empty
        if(!g->vis[id]){
          g->vis[id]=1;
          if(top==cap){ cap*=2; stack=(Node*)realloc(stack,cap*sizeof(Node)); }
          stack[top++] = (Node){nx,ny,nz};
        }
      }
    }
  }
  free(stack);
  return count;
}

// ---------- Helpers ----------
#ifdef __linux__
static double read_rss_mb(){
  struct rusage ru; getrusage(RUSAGE_SELF, &ru);
  // ru_maxrss: kilobytes on Linux
  return ru.ru_maxrss / 1024.0;
}
static void print_cpu_model(){
  FILE* f=fopen("/proc/cpuinfo","r");
  if(!f){ printf("# CPU: unknown\n"); return; }
  char line[512];
  while(fgets(line,sizeof(line),f)){
    if(strncmp(line,"model name",10)==0){
      char *p=strchr(line,':'); if(p){ p++; while(*p==' '){p++;} }
      printf("# CPU: %s", p?p:"unknown\n");
      fclose(f); return;
    }
  }
  fclose(f);
  printf("# CPU: unknown\n");
}
#else
static double read_rss_mb(){ return 0.0; }
static void print_cpu_model(){ printf("# CPU: unknown\n"); }
#endif

typedef struct {
  long long seeds, hits, refined;
  double secs;
} RunStats;

// Core run for one grid
static RunStats run_bls_once(const Config *cfg, int NX, int NY, int NZ){
  // set RNG per run to make comparable
  uint64_t rng = cfg->seed;
  Grid G; grid_init(&G, NX,NY,NZ);
  grid_fill_random(&G, cfg->prob, &rng);

  // basis
  Mat3 B0;
  if(!strcmp(cfg->lattice,"cubic"))              B0 = basis_cubic_unit();
  else if(!strcmp(cfg->lattice,"hexagonal"))     B0 = basis_hex_unit(cfg->hex_c_over_a);
  else if(!strcmp(cfg->lattice,"triclinic"))     B0 = basis_triclinic_unit(cfg->ta,cfg->tb,cfg->tc,
                                                                           cfg->talpha_deg,cfg->tbeta_deg,cfg->tgamma_deg);
  else { fprintf(stderr,"Unknown lattice.\n"); exit(1); }

  Offsets O = centering_offsets(cfg->centering);

  // dNN
  double dNN = cfg->dNN;
  if(dNN<=0.0){
    if(cfg->n_radii<=0){ fprintf(stderr,"Auto dNN but no --radii given. Provide --dnn or --radii.\n"); exit(1); }
    double min_sum=1e18;
    for(int i=0;i<cfg->n_radii;i++)
      for(int j=0;j<cfg->n_radii;j++){
        double s=cfg->radii[i]+cfg->radii[j];
        if(s<min_sum) min_sum=s;
      }
    dNN = cfg->alpha * min_sum;
  }

  double dmin_unit = shortest_unit_spacing(B0,O);
  if(dmin_unit<=1e-12){ fprintf(stderr,"Degenerate unit basis.\n"); exit(1); }
  double scale = dNN / dmin_unit;
  Mat3 B = scale_basis(B0, scale);

  long long seeds_total=0, hits=0, refined_voxels=0;

  clock_t t0=clock();

  for(int oi=0; oi<O.n; ++oi){
    Vec3 ofs = O.ofs[oi];
    Vec3 b1 = M_col(&B,0), b2=M_col(&B,1), b3=M_col(&B,2);
    double bx = fabs(b1.v[0]) + fabs(b2.v[0]) + fabs(b3.v[0]) + 1e-9;
    double by = fabs(b1.v[1]) + fabs(b2.v[1]) + fabs(b3.v[1]) + 1e-9;
    double bz = fabs(b1.v[2]) + fabs(b2.v[2]) + fabs(b3.v[2]) + 1e-9;
    int Imax = (int)ceil((NX + 2.0) / bx) + 2;
    int Jmax = (int)ceil((NY + 2.0) / by) + 2;
    int Kmax = (int)ceil((NZ + 2.0) / bz) + 2;
    if(Imax<0) Imax=NX;
    if(Jmax<0) Jmax=NY;
    if(Kmax<0) Kmax=NZ;

    for(int i=0;i<=Imax;++i){
      for(int j=0;j<=Jmax;++j){
        for(int k=0;k<=Kmax;++k){
          Vec3 p = v3add( M_mul(&B,i,j,k),
                          v3( ofs.v[0]*b1.v[0] + ofs.v[1]*b2.v[0] + ofs.v[2]*b3.v[0],
                              ofs.v[0]*b1.v[1] + ofs.v[1]*b2.v[1] + ofs.v[2]*b3.v[1],
                              ofs.v[0]*b1.v[2] + ofs.v[1]*b2.v[2] + ofs.v[2]*b3.v[2] ) );
          int x = (int)llround(p.v[0]);
          int y = (int)llround(p.v[1]);
          int z = (int)llround(p.v[2]);
          if((unsigned)x >= (unsigned)NX || (unsigned)y >= (unsigned)NY || (unsigned)z >= (unsigned)NZ) continue;
          seeds_total++;
          size_t id = idx3(&G,x,y,z);
          if(G.occ[id] && !G.vis[id]){
            hits++;
            refined_voxels += refine_cluster(&G, x,y,z, cfg->skip);
          }
        }
      }
    }
  }

  clock_t t1=clock();
  free(O.ofs);
  double secs = (double)(t1-t0)/CLOCKS_PER_SEC;

  grid_free(&G);
  RunStats rs = {seeds_total, hits, refined_voxels, secs};
  return rs;
}

// ---------- Main ----------
static double alloc_bytes_MB(int NX,int NY,int NZ){
  double bytes = 2.0 * (double)NX * (double)NY * (double)NZ; // occ + vis, each 1 byte
  return bytes / (1024.0*1024.0);
}

int main(int argc, char**argv){
  Config cfg; parse_cli(argc,argv,&cfg);

  print_cpu_model();
  printf("# lattice=%s centering=%s prob=%.6f seed=%llu alpha=%.3f skip=%d\n",
         cfg.lattice, cfg.centering, cfg.prob, (unsigned long long)cfg.seed, cfg.alpha, cfg.skip);
  if(cfg.dNN>0.0) printf("# dNN=%.6f (explicit)\n", cfg.dNN);
  else            printf("# dNN=alpha*min(r_i+r_j)  (auto)\n");

  if(cfg.bench){
    printf("N\tvoxels\talloc_MB\tmaxRSS_MB\ttime_s\tseeds\thits\n");
    for(int N=cfg.nmin; N<=cfg.nmax; N+=cfg.nstep){
      // square box
      RunStats rs = run_bls_once(&cfg, N, N, N);
      double rssMB = read_rss_mb();
      double allocMB = alloc_bytes_MB(N,N,N);
      long long vox = (long long)N*(long long)N*(long long)N;
      printf("%d\t%lld\t%.3f\t%.3f\t%.6f\t%lld\t%lld\n",
             N, vox, allocMB, rssMB, rs.secs, rs.seeds, rs.hits);
      fflush(stdout);
    }
    return 0;
  }

  // single-run
  RunStats rs = run_bls_once(&cfg, cfg.NX, cfg.NY, cfg.NZ);
  double rssMB = read_rss_mb();
  double allocMB = alloc_bytes_MB(cfg.NX,cfg.NY,cfg.NZ);

  printf("# single run\n");
  printf("grid=%dx%dx%d prob=%.6f\n", cfg.NX,cfg.NY,cfg.NZ, cfg.prob);
  printf("alloc_MB=%.3f maxRSS_MB=%.3f time_s=%.6f seeds=%lld hits=%lld refined_voxels=%lld\n",
         allocMB, rssMB, rs.secs, rs.seeds, rs.hits, rs.refined);

  if(cfg.radii) free(cfg.radii);
  return 0;
}

