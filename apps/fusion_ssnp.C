#define VERSIONED 1
#define WEIGHTED 1
#include "ligra.h"

using IdxType = long long;
struct DJ_F {
  intE* NarrowestPathVal;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;
  intE versionNum;

  DJ_F(intE* _NarrowestPathVal, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize, intE _versionNum) : 
    NarrowestPathVal(_NarrowestPathVal), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize), versionNum(_versionNum) {}
  
  inline bool update(uintE s, uintE d, intE edgeLen, intE versionTag) 
  { 
    bool ret = false;

    for (int v = 0; v < versionNum; v++) 
    {  
      if (!graphUtils::validate(v, (versionNum - 1), versionTag)) 
      {
        continue;
      }
      for (long j = 0; j < BatchSize; j++) 
      {
        IdxType s_index = s * (versionNum * BatchSize) + v + versionNum * j;
        IdxType d_index = d * (versionNum * BatchSize) + v + versionNum * j;

        if (CurrActiveArray[s_index]) {
          intE newValue = std::max(NarrowestPathVal[s_index], edgeLen);
          if (NarrowestPathVal[d_index] > newValue) 
          {
            NarrowestPathVal[d_index] = newValue;
            NextActiveArray[d_index] = true;
            ret = true;
          }
        }
      }
    }
    return ret;
  }

  inline bool updateAtomic(uintE s, uintE d, intE edgeLen, intE versionTag) { 
    bool ret = false;

    for (long j = 0; j < BatchSize; j++) {  
      for (int v = 0; v < versionNum; v++) {  
        if (!graphUtils::validate(v, (versionNum - 1), versionTag)) {
          continue;
        }
        IdxType s_index = s * (versionNum * BatchSize) + v + versionNum * j;
        IdxType d_index = d * (versionNum * BatchSize) + v + versionNum * j;

        if (CurrActiveArray[s_index]) 
        {
          intE newValue = std::max(NarrowestPathVal[s_index], edgeLen);
          
          if (writeMin(&NarrowestPathVal[d_index], newValue)) {
            if (CAS(&NextActiveArray[d_index], false, true)) {
              ret = true;
            }
          }
        }
      }
    }
    return ret;
  }

  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SINGLE_F {
  intE* NarrowestPathVal;
  intE targetVersion;
  intE versionNum;
  DJ_SINGLE_F(intE* _NarrowestPathVal, intE _targetVersion, intE _versionNum) : 
    NarrowestPathVal(_NarrowestPathVal), targetVersion(_targetVersion), versionNum(_versionNum) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag) { 
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }    
    intE newValue = std::max(NarrowestPathVal[s], edgeLen);
    if (NarrowestPathVal[d] > newValue) {
      NarrowestPathVal[d] = newValue;
      ret = true;
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag){ 
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }
    IdxType s_begin = s;
    IdxType d_begin = d;
    intE newValue = std::max(NarrowestPathVal[s], edgeLen);
    if (writeMin(&NarrowestPathVal[d], newValue)) {
      ret = true;
    }
    return ret;
  }
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SKIP_F {
  intE* NarrowestPathVal;
  long BatchSize;
  intE totalVersion;
  intE parallelVersionNum;
  intE* versionNum;

  DJ_SKIP_F(intE* _NarrowestPathVal, 
            long _BatchSize, 
            intE _totalVersion = 0,
            intE _parallelVersionNum = 0, 
            intE* _versionNum = nullptr) : 
    NarrowestPathVal(_NarrowestPathVal), BatchSize(_BatchSize),
    totalVersion(_totalVersion),  
    parallelVersionNum(_parallelVersionNum), versionNum(_versionNum) {}
  
  inline bool update(uintE s, uintE d, intE edgeLen, intE versionTag = commonTag) 
  {
    bool ret = false;
  
    for (long j = 0; j < BatchSize; j++) {
      for (intE v = 0; v < parallelVersionNum; v++) {
        intE real_version = versionNum[v];
        if (!graphUtils::validate(real_version, totalVersion - 1, versionTag)) {
          continue;
        }
        IdxType s_idx = ((IdxType)s * BatchSize + j) * parallelVersionNum + v;
        IdxType d_idx = ((IdxType)d * BatchSize + j) * parallelVersionNum + v;
        intE newValue = std::max(NarrowestPathVal[s_idx], edgeLen);
        if (NarrowestPathVal[d_idx] > newValue) {
          NarrowestPathVal[d_idx] = newValue;
          ret = true;
        }
      }
    }
    return ret;
  }

  inline bool updateAtomic(uintE s, uintE d, intE edgeLen, intE versionTag = commonTag) {
    bool ret = false;
    for (long j = 0; j < BatchSize; j++) {
      for (intE v = 0; v < parallelVersionNum; v++) {
        intE real_version = versionNum[v];
        if (!graphUtils::validate(real_version, totalVersion - 1, versionTag)) 
        {
          continue;
        }
        IdxType s_idx = ((IdxType)s * BatchSize + j) * parallelVersionNum + v;
        IdxType d_idx = ((IdxType)d * BatchSize + j) * parallelVersionNum + v;

        intE newValue = std::max(NarrowestPathVal[s_idx], edgeLen);

        if (writeMin(&NarrowestPathVal[d_idx], newValue)) {
          ret = true;
        }
      }
    }
    return ret;
  }

  inline bool cond (uintE d) { return cond_true(d); } 
};

struct BFSLV_F {
  uintE* Levels;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  BFSLV_F(uintE* _Levels, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    Levels(_Levels), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen=1) { return false; }
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen=1) { return false; }
  inline bool cond (uintE d) { return false; } 
};

template <class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  intE numVersions = G.batchNum + 1;
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)numVersions;
  fprintf(stderr, "totalNumVertices: %lld\n", totalNumVertices);
  intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for (long i = 0; i < batch_size; i++) 
  { 
    for (int v = 0; v < numVersions; v++) 
    {
        int idx = vecQueries[i] * (numVersions * batch_size) + v + numVersions * i;
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
    }
  }

  for (long i = 0; i < batch_size; i++) 
  {  
    for (int v = 0; v < numVersions; v++) 
    {
        int idx = vecQueries[i] * (numVersions * batch_size) + v + numVersions * i;
        CurrActiveArray[idx] = true;
    }
  }

  vertexSubset Frontier(n, frontier);

  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    vertexSubset output = edgeMap(G, Frontier, DJ_F(NarrowestPathVal, CurrActiveArray, NextActiveArray, batch_size, numVersions), -1, no_dense|remove_duplicates);
    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }
  cout << "Total iterations: " << iteration << endl;

  Frontier.del();
  pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>& G, std::vector<long> vecQueries, queue<long>& queryQueue, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    NarrowestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;  // Set source nodes to 0
  }
  for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);

  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  int slots_parameter = P.getOptionIntValue("-slots", 1);
  int empty_slots = 0;
  vector<long> next_batch;
  vector<long> lastQueries;
  bool _flag_check = true;
  bool* finish_status = pbbs::new_array<bool>(batch_size);
  for(long i = 0; i < batch_size; i++) {
    finish_status[i] = true;
  }
  
  vector<long> queuedQueries;
  long queued_index = 0;
  while (!queryQueue.empty()) {
    queuedQueries.push_back(queryQueue.front());
    queryQueue.pop();
  }

  while(!Frontier.isEmpty() || queued_index != queuedQueries.size()){
    iteration++;
    totalActivated += Frontier.size();

    vertexSubset output = edgeMap(G, Frontier, DJ_F(NarrowestPathVal, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }

    bool* is_batch_finished = pbbs::new_array<bool>(batch_size);
    for(size_t i = 0; i < batch_size; i++) {
      is_batch_finished[i] = false;
    }
    empty_slots = 0;
    if (queued_index != queuedQueries.size()) {
      vector<long> empty_slots_vec;
      for (int i = 0; i < batch_size; i++) {
        bool is_finished = true;
        long active_cnt = 0;
        parallel_for(size_t index = 0; index < n; index++) {
          if (CurrActiveArray[index * batch_size + i]) {
            pbbs::fetch_and_add(&active_cnt, 1);
          }
        }
        if (active_cnt == 0) {
          empty_slots++;
          empty_slots_vec.push_back(i);
          if (finish_status[i]) {
            finish_status[i] = false;
          }
        }
      }
      if (empty_slots == batch_size) {
        cout << "batch is empty, iteration " << iteration << endl;
      }

      if (empty_slots >= slots_parameter || (queuedQueries.size()-queued_index <= empty_slots)) {
        for (int i = 0; i < empty_slots_vec.size(); i++) {
          if (queued_index < queuedQueries.size()) {
            long nextQ = queuedQueries[queued_index];
            queued_index++;
            int index_in_batch = empty_slots_vec[i];
            vecQueries[index_in_batch] = nextQ;
            finish_status[index_in_batch] = true;

            CurrActiveArray[nextQ * batch_size + index_in_batch] = true;
            parallel_for(size_t index = 0; index < n; index++) {
              NarrowestPathVal[index * batch_size + index_in_batch] = (intE)MAXPATH;  // Initialize to MAXPATH
            }
            NarrowestPathVal[nextQ * batch_size + index_in_batch] = 0;  // Set source to 0

            Frontier.toDense();
            bool* new_d = pbbs::new_array<bool>(n);
            parallel_for(size_t ii = 0; ii < n; ii++) {
              new_d[ii] = Frontier.d[ii];
            }
            new_d[nextQ] = true;
            vertexSubset Frontier_new(n, new_d);
            Frontier.del();
            Frontier = Frontier_new;
          }
        }
      }
    }
  }

  Frontier.del();
  pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Skipping(graph<vertex>& G, 
  std::vector<long> vecQueries, commandLine P, bool should_profile) 
{
#ifdef VERSIONED
  auto totalVersion = G.batchNum;
  std::cout << "Number of Nodes in the graph is " << G.n << std::endl;
  std::cout << "Number of Edges in the graph is " << G.m << std::endl;
  std::cout << "Total Number of Version is " << totalVersion << std::endl;
  auto parallelVersionNum = P.getOptionIntValue("-parallelVersion", 4);
  double total_propagation_time = 0.0;
  double baseTime = 0.0;

  for (intE versionStart = 0; versionStart < totalVersion; versionStart += parallelVersionNum) {
    intE versionNum = std::min<intE>(parallelVersionNum, totalVersion - versionStart);
    intE* versionNumArray = newA(intE, versionNum);
    for (intE i = 0; i < versionNum; i++) {
      versionNumArray[i] = versionStart + i;
    }
#endif

    size_t n = G.n;
    size_t edge_count = G.m;
    long batch_size = vecQueries.size();
#ifndef VERSIONED
    IdxType totalNumVertices = (IdxType)n * (IdxType)versionNum;
#else
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;
#endif
    intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    }

#ifndef VERSIONED
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * versionNum + v);
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }
#else
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }
#endif

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;
    size_t totalNoOverlap = 0;

#ifdef VERSIONED
    timer t_propagate;
    t_propagate.start();
#endif

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();
#ifndef VERSIONED
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(NarrowestPathVal, batch_size), -1, no_dense | remove_duplicates);
#else
      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(NarrowestPathVal, batch_size, totalVersion, versionNum, versionNumArray),
        -1, no_dense | remove_duplicates);
#endif
      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;
      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }

#ifdef VERSIONED
    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;
#endif
    Frontier.del();

#ifdef VERSIONED
    for (intE query = 0; query < batch_size; query++) {
      intE src = vecQueries[query];
      for (intE v = 0; v < versionNum; v++) {
        intE version = versionNumArray[v];
        intE* tmpVal = pbbs::new_array<intE>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = (intE)MAXPATH;  // Initialize to MAXPATH
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;  // Set source to 0
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            DJ_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  }
#ifdef VERSIONED
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
  std::cout << "Total baseline time: " << baseTime << " seconds\n";
  std::cout << "Speedup: " << baseTime / total_propagation_time << std::endl;
#endif
  return make_pair(0, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      frontier[vecQueries[i]] = true;
    }
  }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    NarrowestPathVal[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;  // Set source nodes to 0
  }
  for(size_t i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
    }
  }

  vertexSubset Frontier(n, frontier);

  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    vertexSubset output = edgeMap(G, Frontier, DJ_F(NarrowestPathVal, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
        NextActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
      }
    }
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G,
  std::vector<long> vecQueries, commandLine P,
  std::vector<int> defer_vec, bool should_profile)
{
#ifdef VERSIONED
  auto totalVersion = G.batchNum;
  auto parallelVersionNum = P.getOptionIntValue("-parallelVersion", 4);
  double total_propagation_time = 0.0;
  double baseTime = 0.0;

  for (intE versionStart = 0; versionStart < totalVersion; versionStart += parallelVersionNum) {
    intE versionNum = std::min<intE>(parallelVersionNum, totalVersion - versionStart);
    intE* versionNumArray = newA(intE, versionNum);
    for (intE i = 0; i < versionNum; i++) {
      versionNumArray[i] = versionStart + i;
    }
#endif

    size_t n = G.n;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) {
      frontier[i] = false;
    }

    for (long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == 0) {
        frontier[vecQueries[i]] = true;
      }
    }

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;

#ifdef VERSIONED
    timer t_propagate;
    t_propagate.start();
#endif

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();

      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(NarrowestPathVal, batch_size, totalVersion, versionNum, versionNumArray),
        -1, no_dense | remove_duplicates);

      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;

      for (long i = 0; i < batch_size; i++) {
        if (defer_vec[i] == iteration) {
          new_d[vecQueries[i]] = true;
        }
      }

      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }

#ifdef VERSIONED
    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;
#endif

    Frontier.del();
#ifdef VERSIONED
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  }
#ifdef VERSIONED
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
#endif
  return make_pair(0, 0);
}

template <class vertex>
pair<double, double> Compute_Base_Skipping_Time(graph<vertex>& G, 
  std::vector<long> vecQueries, commandLine P, bool should_profile) 
{
#ifdef VERSIONED
  auto totalVersion = G.batchNum;
  std::cout << "Number of Nodes in the graph is " << G.n << std::endl;
  std::cout << "Number of Edges in the graph is " << G.m << std::endl;
  std::cout << "Total Number of Version is " << totalVersion << std::endl;
  auto parallelVersionNum = P.getOptionIntValue("-parallelVersion", 4);
  double total_propagation_time = 0.0;
  double baseTime = 0.0;

  for (intE versionStart = 0; versionStart < totalVersion; versionStart += parallelVersionNum) {
    intE versionNum = std::min<intE>(parallelVersionNum, totalVersion - versionStart);
    intE* versionNumArray = newA(intE, versionNum);
    for (intE i = 0; i < versionNum; i++) {
      versionNumArray[i] = versionStart + i;
    }
#endif

    size_t n = G.n;
    size_t edge_count = G.m;
    long batch_size = vecQueries.size();
#ifndef VERSIONED
    IdxType totalNumVertices = (IdxType)n * (IdxType)versionNum;
#else
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;
#endif
    intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    }

#ifndef VERSIONED
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * versionNum + v);
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }
#else
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }
#endif

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;
    size_t totalNoOverlap = 0;

#ifdef VERSIONED
    timer t_propagate;
    t_propagate.start();
#endif

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();
#ifndef VERSIONED
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(NarrowestPathVal, batch_size), -1, no_dense | remove_duplicates);
#else
      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(NarrowestPathVal, batch_size, totalVersion, versionNum, versionNumArray),
        -1, no_dense | remove_duplicates);
#endif
      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;
      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }

#ifdef VERSIONED
    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;
#endif
    Frontier.del();

#ifdef VERSIONED
    for (intE query = 0; query < batch_size; query++) {
      intE src = vecQueries[query];
      for (intE v = 0; v < versionNum; v++) {
        intE version = versionNumArray[v];
        intE* tmpVal = pbbs::new_array<intE>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = (intE)MAXPATH;  // Initialize to MAXPATH
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;  // Set source to 0
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            DJ_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * batch_size + query) * versionNum + v;
        //   if (NarrowestPathVal[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << version << " error: " << error << std::endl;          
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  }
#ifdef VERSIONED
  return make_pair(total_propagation_time, baseTime);
#else
  return make_pair(0.0, 0.0);
#endif  
}

template <class vertex>
pair<double, double> Compute_Delay_Skipping_Time(graph<vertex>& G,
  std::vector<long> vecQueries, commandLine P,
  std::vector<int> defer_vec, bool should_profile)
{
#ifdef VERSIONED
  auto totalVersion = G.batchNum;
  auto parallelVersionNum = P.getOptionIntValue("-parallelVersion", 4);
  double total_propagation_time = 0.0;
  double baseTime = 0.0;

  for (intE versionStart = 0; versionStart < totalVersion; versionStart += parallelVersionNum) {
    intE versionNum = std::min<intE>(parallelVersionNum, totalVersion - versionStart);
    intE* versionNumArray = newA(intE, versionNum);
    for (intE i = 0; i < versionNum; i++) {
      versionNumArray[i] = versionStart + i;
    }
#endif

    size_t n = G.n;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    intE* NarrowestPathVal = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) {
      frontier[i] = false;
    }

    for (long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == 0) {
        frontier[vecQueries[i]] = true;
      }
    }

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NarrowestPathVal[i] = (intE)MAXPATH;  // Initialize all nodes to MAXPATH
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        NarrowestPathVal[idx] = 0;  // Set source nodes to 0
      }
    }

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;

#ifdef VERSIONED
    timer t_propagate;
    t_propagate.start();
#endif

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();

      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(NarrowestPathVal, batch_size, totalVersion, versionNum, versionNumArray),
        -1, no_dense | remove_duplicates);

      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;

      for (long i = 0; i < batch_size; i++) {
        if (defer_vec[i] == iteration) {
          new_d[vecQueries[i]] = true;
        }
      }

      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }

#ifdef VERSIONED
    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;
#endif

    Frontier.del();

#ifdef VERSIONED
    for (intE query = 0; query < batch_size; query++) {
      intE src = vecQueries[query];
      for (intE v = 0; v < versionNum; v++) {
        intE version = versionNumArray[v];
        intE* tmpVal = pbbs::new_array<intE>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = (intE)MAXPATH;  // Initialize to MAXPATH
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;  // Set source to 0
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            DJ_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * batch_size + query) * versionNum + v;
        //   if (NarrowestPathVal[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << version << " error: " << error << std::endl;          
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(NarrowestPathVal, totalNumVertices);
  }
#ifdef VERSIONED
  return make_pair(total_propagation_time, baseTime);
#else
  return make_pair(0.0, 0.0);
#endif  
}

template <class vertex>
uintE* Compute_Eval(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Hops = pbbs::new_array<uintE>(totalNumVertices);
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Hops[i] = 0;  // 初始化为0
  }
  return Hops;
} 