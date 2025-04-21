#define VERSIONED 1
#define WEIGHTED 1
#include "ligra.h"

using IdxType = long long;

struct SSR_F {
  uintE* Reachable;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;
  intE versionNum;

  SSR_F(uintE* _Reachable, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize, intE _versionNum) : 
    Reachable(_Reachable), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize), versionNum(_versionNum) {}
  
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
        IdxType s_index = ((IdxType)s * BatchSize + j) * versionNum + v;
        IdxType d_index = ((IdxType)d * BatchSize + j) * versionNum + v;

        if (CurrActiveArray[s_index]) {
          uintE newValue = (Reachable[s_index] > 0 || Reachable[d_index] > 0) ? 1 : 0;
          if (Reachable[d_index] != newValue) 
          {
            Reachable[d_index] = newValue;
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
        IdxType s_index = ((IdxType)s * BatchSize + j) * versionNum + v;
        IdxType d_index = ((IdxType)d * BatchSize + j) * versionNum + v;

        if (CurrActiveArray[s_index]) 
        {
          uintE newValue = (Reachable[s_index] > 0 || Reachable[d_index] > 0) ? 1 : 0;
          uintE oldValue = (newValue == 0) ? 1 : 0;
          if (CAS(&Reachable[d_index], oldValue, newValue)) {
            if (CAS(&NextActiveArray[d_index], false, true)) {
              ret = true;
            }
          }
        }
      }
    }
    return ret;
  }

  inline bool cond(uintE d) { return cond_true(d); }
};

struct SSR_SINGLE_F {
  bool* Reachable;
  intE targetVersion;
  intE versionNum;
  SSR_SINGLE_F(bool* _Reachable, intE _targetVersion, intE _versionNum) : 
    Reachable(_Reachable), targetVersion(_targetVersion), versionNum(_versionNum) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag) { 
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }    
    bool newValue = Reachable[s] || Reachable[d];
    if (Reachable[d] != newValue) {
      Reachable[d] = newValue;
      ret = true;
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag){ 
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }
    bool newValue = Reachable[s] || Reachable[d];
    if (CAS(&Reachable[d], !newValue, newValue)) {
      ret = true;
    }
    return ret;
  }
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct SSR_SKIP_F {
  bool* Reachable;
  long BatchSize;
  intE totalVersion;
  intE parallelVersionNum;
  intE* versionNum;

  SSR_SKIP_F(bool* _Reachable, 
            long _BatchSize, 
            intE _totalVersion = 0,
            intE _parallelVersionNum = 0, 
            intE* _versionNum = nullptr) : 
    Reachable(_Reachable), BatchSize(_BatchSize),
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
        bool newValue = Reachable[s_idx] || Reachable[d_idx];
        if (Reachable[d_idx] != newValue) {
          Reachable[d_idx] = newValue;
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

        bool newValue = Reachable[s_idx] || Reachable[d_idx];

        if (CAS(&Reachable[d_idx], !newValue, newValue)) {
          ret = true;
        }
      }
    }
    return ret;
  }

  inline bool cond (uintE d) { return cond_true(d); } 
};

template <class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  intE numVersions = G.batchNum + 1;
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)numVersions;
  fprintf(stderr, "totalNumVertices: %lld\n", totalNumVertices);
  uintE* Reachable = pbbs::new_array<uintE>(totalNumVertices);
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
    Reachable[i] = 0;  // Initialize all nodes to 0
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for (long i = 0; i < batch_size; i++) 
  { 
    for (int v = 0; v < numVersions; v++) 
    {
        int idx = vecQueries[i] * (numVersions * batch_size) + v + numVersions * i;
        Reachable[idx] = 1;  // Set source nodes to 1
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

    vertexSubset output = edgeMap(G, Frontier, SSR_F(Reachable, CurrActiveArray, NextActiveArray, batch_size, numVersions), -1, no_dense|remove_duplicates);
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
  pbbs::delete_array(Reachable, totalNumVertices);
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
    bool* Reachable = pbbs::new_array<bool>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      Reachable[i] = false;  // Initialize all nodes to false
    }

#ifndef VERSIONED
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * versionNum + v);
        Reachable[idx] = true;  // Set source nodes to true
      }
    }
#else
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        Reachable[idx] = true;  // Set source nodes to true
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
      vertexSubset output = edgeMap(G, Frontier, SSR_SKIP_F(Reachable, batch_size), -1, no_dense | remove_duplicates);
#else
      vertexSubset output = edgeMap(G, Frontier,
        SSR_SKIP_F(Reachable, batch_size, totalVersion, versionNum, versionNumArray),
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
        bool* tmpVal = pbbs::new_array<bool>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = false;  // Initialize to false
          tmpFrontier[i] = false;
        }

        tmpVal[src] = true;  // Set source to true
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            SSR_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * batch_size * versionNum) + 
        //                ((IdxType)query * versionNum) + v;
        //   if (Reachable[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << versionNumArray[v] << " error: " << error << std::endl;          
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(Reachable, totalNumVertices);
  }
#ifdef VERSIONED
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
  std::cout << "Total baseline time: " << baseTime << " seconds\n";
  std::cout << "Speedup: " << baseTime / total_propagation_time << std::endl;
#endif
  return make_pair(0, 0);
}

template <class vertex>
uintE* Compute_Eval(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Reachable = pbbs::new_array<uintE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  
  // Initialize frontier
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    frontier[vecQueries[i]] = true;
  }

  // Initialize arrays
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Reachable[i] = 0;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }

  // Set source nodes to 1
  for(long i = 0; i < batch_size; i++) {
    Reachable[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 1;
  }

  // Set source nodes as active
  parallel_for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);
  while(!Frontier.isEmpty()){
    vertexSubset output = edgeMap(G, Frontier, SSR_F(Reachable, CurrActiveArray, NextActiveArray, batch_size, 1), -1, no_dense|remove_duplicates);
    Frontier.del();
    Frontier = output;

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return Reachable;
}

template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G,
  std::vector<long> vecQueries, commandLine P,
  std::vector<int> defer_vec, bool should_profile)
{
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  bool* Reachable = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0 || defer_vec[i] > 1000) {
      frontier[vecQueries[i]] = true;
    }
  }

  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    Reachable[i] = false;  // Initialize all nodes to false
  }
  for(long i = 0; i < batch_size; i++) {
    Reachable[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;  // Set source nodes to true
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    vertexSubset output = edgeMap(G, Frontier, SSR_SKIP_F(Reachable, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    Frontier.toDense();
    bool* new_d = Frontier.d;
    Frontier.d = nullptr;
    for(long i = 0; i < batch_size; i++) {
      if (defer_vec[i] == iteration) {
        new_d[vecQueries[i]] = true;
      } 
    }
    vertexSubset Frontier_new(n, new_d);
    Frontier.del();
    Frontier = Frontier_new;
  }

  Frontier.del();
  pbbs::delete_array(Reachable, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<double, double> Compute_Delay_Skipping_Time(graph<vertex>& G,
  std::vector<long> vecQueries, commandLine P,
  std::vector<int> defer_vec, bool should_profile)
{
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

    size_t n = G.n;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    bool* Reachable = pbbs::new_array<bool>(totalNumVertices);
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
      Reachable[i] = false;  // Initialize all nodes to false
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * (IdxType)batch_size * (IdxType)versionNum) + 
                     ((IdxType)i * (IdxType)versionNum) + v;
        Reachable[idx] = true;  // Set source nodes to true
      }
    }

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;

    timer t_propagate;
    t_propagate.start();

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();

      vertexSubset output = edgeMap(G, Frontier,
        SSR_SKIP_F(Reachable, batch_size, totalVersion, versionNum, versionNumArray),
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

    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;

    Frontier.del();

    for (intE query = 0; query < batch_size; query++) {
      intE src = vecQueries[query];
      for (intE v = 0; v < versionNum; v++) {
        intE version = versionNumArray[v];
        bool* tmpVal = pbbs::new_array<bool>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = false;  // Initialize to false
          tmpFrontier[i] = false;
        }

        tmpVal[src] = true;  // Set source to true
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            SSR_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * (IdxType)batch_size * (IdxType)versionNum) + 
        //                ((IdxType)query * (IdxType)versionNum) + v;
        //   if (Reachable[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << versionNumArray[v] << " error: " << error << std::endl;          
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
    pbbs::delete_array(Reachable, totalNumVertices);
  }
  return make_pair(total_propagation_time, baseTime);
}

template <class vertex>
pair<double, double> Compute_Base_Skipping_Time(graph<vertex>& G, 
  std::vector<long> vecQueries, commandLine P, bool should_profile) 
{
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

    size_t n = G.n;
    size_t edge_count = G.m;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    bool* Reachable = pbbs::new_array<bool>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      Reachable[i] = false;  // Initialize all nodes to false
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * (IdxType)batch_size * (IdxType)versionNum) + 
                     ((IdxType)i * (IdxType)versionNum) + v;
        Reachable[idx] = true;  // Set source nodes to true
      }
    }

    vertexSubset Frontier(n, frontier);

    long iteration = 0;
    size_t totalActivated = 0;
    size_t totalNoOverlap = 0;

    timer t_propagate;
    t_propagate.start();

    while (!Frontier.isEmpty()) {
      iteration++;
      totalActivated += Frontier.size();

      vertexSubset output = edgeMap(G, Frontier,
        SSR_SKIP_F(Reachable, batch_size, totalVersion, versionNum, versionNumArray),
        -1, no_dense | remove_duplicates);

      Frontier.del();
      Frontier = output;

      Frontier.toDense();
      bool* new_d = Frontier.d;
      Frontier.d = nullptr;
      vertexSubset Frontier_new(n, new_d);
      Frontier.del();
      Frontier = Frontier_new;
    }

    t_propagate.stop();
    total_propagation_time += t_propagate.totalTime;

    Frontier.del();

    for (intE query = 0; query < batch_size; query++) {
      intE src = vecQueries[query];
      for (intE v = 0; v < versionNum; v++) {
        intE version = versionNumArray[v];
        bool* tmpVal = pbbs::new_array<bool>(n);
        bool* tmpFrontier = pbbs::new_array<bool>(n);

        parallel_for(size_t i = 0; i < n; i++) {
          tmpVal[i] = false;  // Initialize to false
          tmpFrontier[i] = false;
        }

        tmpVal[src] = true;  // Set source to true
        tmpFrontier[src] = true;
        vertexSubset Frontier_tmp(n, tmpFrontier);
        timer baseTimer;
        baseTimer.start();
        while (!Frontier_tmp.isEmpty()) {
          vertexSubset output = edgeMap(G, Frontier_tmp,
            SSR_SINGLE_F(tmpVal, version, totalVersion),
            -1, no_dense | remove_duplicates);
          Frontier_tmp.del();
          Frontier_tmp = output;
        }
        baseTimer.stop();
        baseTime += baseTimer.totalTime;
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * (IdxType)batch_size * (IdxType)versionNum) + 
        //                ((IdxType)query * (IdxType)versionNum) + v;
        //   if (Reachable[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << versionNumArray[v] << " error: " << error << std::endl;          
        Frontier_tmp.del();
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
    pbbs::delete_array(Reachable, totalNumVertices);
  }
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
  std::cout << "Total baseline time: " << baseTime << " seconds\n";
  std::cout << "Speedup: " << baseTime / total_propagation_time << std::endl;
  return make_pair(total_propagation_time, baseTime);
} 