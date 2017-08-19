#include "../../include/smctc.hh"

#define MAXVAR 20
#define NUMNODEVAR 2

union INDEX_TYPE_t {
  unsigned short var[MAXVAR];
  unsigned long long index;
};

typedef INDEX_TYPE_t INDEX_TYPE;

class cv_state
{
public:
  INDEX_TYPE index;
  int size;
  double logBF;
};

double get_logBF(INDEX_TYPE index,int size);

inline bool is_in_model(const INDEX_TYPE& I,int i);
inline bool is_not_in_model(const INDEX_TYPE& I, int i);

void get_model_data(INDEX_TYPE index, int size);
double compute_rho(int level, int rho_method = 0);



smc::particle<cv_state> fInitialise(smc::rng *pRng);

long fSelect(long lTime, const smc::particle<cv_state> & p,
	     smc::rng *pRng);

void fMove(long lTime, smc::particle<cv_state> & pFrom,
	   smc::rng *pRng);

void fMove2(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);
void fMove3(long lTime, smc::particle<cv_state > & pFrom, smc::rng *pRng);


double get_add_prob(cv_state *model_ptr, int i);
inline int get_i_sub(INDEX_TYPE I, int i);
void build_subGBT(const cv_state & root_model);
INDEX_TYPE make_global_index(INDEX_TYPE I_root, INDEX_TYPE I_sub, int level);
inline INDEX_TYPE make_child_index(INDEX_TYPE& I, unsigned short add_dim);
inline INDEX_TYPE get_next_node(INDEX_TYPE& I, int n, int k);
double *get_node(INDEX_TYPE& I, int level);
unsigned int get_node_index(INDEX_TYPE& I,int level);
