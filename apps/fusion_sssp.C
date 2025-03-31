// This code is part of the project "Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing"
// Copyright (c) 2022 Xizhe Yin
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// #pragma push_macro("VERSIONED")  // 保存当前的 VERSIONED 状态
// #undef VERSIONED
#define VERSIONED 1  // 仅在 DJ_F 结构体中启用
#define WEIGHTED 1
#include "ligra.h"

using IdxType = long long;
struct DJ_F {
  intE* ShortestPathLen;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;
  intE versionNum;

  // Difference Here is We Need to Update all version result in one update Atomic
  // Array Should BE Like [node0_{version 0}, node0_{version 1}, node0_{version 2}, node1_{version 0}, node1_{version 1}, node1_{version 2}, ...]
  // Each Batch Contains version * batch of query * number of nodes
  DJ_F(intE* _ShortestPathLen, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize, intE _versionNum) : 
    ShortestPathLen(_ShortestPathLen), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize), versionNum(_versionNum) {
      // std::cout << "versionNum: " << versionNum << std::endl;
      // std::cout << "BatchSize: " << BatchSize << std::endl;
    }
  
inline bool update(uintE s, uintE d, intE edgeLen, intE versionTag) 
{ 
  bool ret = false;

  for (int v = 0; v < versionNum; v++) 
  {  // 版本变化最快
    if (!graphUtils::validate(v, (versionNum - 1), versionTag)) 
    {
      continue;
    }
    for (long j = 0; j < BatchSize; j++) 
    {
        IdxType s_index = s * (versionNum * BatchSize) + v + versionNum * j;
        IdxType d_index = d * (versionNum * BatchSize) + v + versionNum * j;

        if (CurrActiveArray[s_index]) {
            intE newValue = ShortestPathLen[s_index] + edgeLen;
            if (ShortestPathLen[d_index] > newValue) 
            {
                ShortestPathLen[d_index] = newValue;
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

  for (long j = 0; j < BatchSize; j++) {  // batch 变化最慢
      for (int v = 0; v < versionNum; v++) {  // version 变化最快
          if (!graphUtils::validate(v, (versionNum - 1), versionTag)) {
              continue;
          }
          IdxType s_index = s * (versionNum * BatchSize) + v + versionNum * j;
          IdxType d_index = d * (versionNum * BatchSize) + v + versionNum * j;

          if (CurrActiveArray[s_index]) 
          {
              intE newValue = ShortestPathLen[s_index] + edgeLen;
              
              if (writeMin(&ShortestPathLen[d_index], newValue)) {
                  if (CAS(&NextActiveArray[d_index], false, true)) {
                      ret = true;
                  }
              }
          }
      }
  }
  return ret;
}

  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SINGLE_F {
  intE* ShortestPathLen;
  intE targetVersion;
  intE versionNum;
  DJ_SINGLE_F(intE* _ShortestPathLen, intE _targetVersion, intE _versionNum) : 
    ShortestPathLen(_ShortestPathLen), targetVersion(_targetVersion), versionNum(_versionNum) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag) { //Update
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }    
    intE newValue = ShortestPathLen[s] + edgeLen;
    if (ShortestPathLen[d] > newValue) {
      ShortestPathLen[d] = newValue;
      ret = true;
    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen, intE versionTag = commonTag){ //atomic version of Update
    bool ret = false;
    if (!graphUtils::validate(targetVersion, versionNum - 1, versionTag)) {
      return ret;
    }
    IdxType s_begin = s;
    IdxType d_begin = d;
    intE newValue = ShortestPathLen[s] + edgeLen;
    if (writeMin(&ShortestPathLen[d], newValue)) {
      ret = true;
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct DJ_SKIP_F {
  intE* ShortestPathLen;
  long BatchSize;
  // The version being computed here
  intE totalVersion; //TotalBatch + 1
  intE parallelVersionNum; // How many versions being batched
  intE* versionNum;

  DJ_SKIP_F(intE* _ShortestPathLen, 
            long _BatchSize, 
            intE _totalVersion = 0,
            intE _parallelVersionNum = 0, 
            intE* _versionNum = nullptr) : 
    ShortestPathLen(_ShortestPathLen), BatchSize(_BatchSize),
    totalVersion(_totalVersion),  
    parallelVersionNum(_parallelVersionNum), versionNum(_versionNum)
    {
// #ifdef VERSIONED
//       std::cout << "versionNum: " << parallelVersionNum << std::endl;
//       std::cout << "totalVersion: " << totalVersion << std::endl;
//       std::cout << "Query BatchSize: " << BatchSize << std::endl;
// #endif
    }
  
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
              intE newValue = ShortestPathLen[s_idx] + edgeLen;
              if (ShortestPathLen[d_idx] > newValue) {
                  ShortestPathLen[d_idx] = newValue;
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

            intE newValue = ShortestPathLen[s_idx] + edgeLen;

            if (writeMin(&ShortestPathLen[d_idx], newValue)) {
                ret = true;
            }
        }
    }
    return ret;
}

  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

struct BFSLV_F {
  uintE* Levels;
  bool* CurrActiveArray;
  bool* NextActiveArray;
  long BatchSize;

  BFSLV_F(uintE* _Levels, bool* _CurrActiveArray, bool* _NextActiveArray, long _BatchSize) : 
    Levels(_Levels), CurrActiveArray(_CurrActiveArray), NextActiveArray(_NextActiveArray), BatchSize(_BatchSize) {}
  
  inline bool update (uintE s, uintE d, intE edgeLen=1, intE versionTag = commonTag) { //Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;

    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        uintE newValue = Levels[s_begin + j] + 1;
        if (Levels[d_begin + j] > newValue) {
          // Visited[d_begin + j] = true;
          Levels[d_begin + j] = newValue;
          NextActiveArray[d_begin + j] = true;
          ret = true;
        }
      }

    }
    return ret;
  }
  
  inline bool updateAtomic (uintE s, uintE d, intE edgeLen=1, intE versionTag = commonTag){ //atomic version of Update
    bool ret = false;
    IdxType s_begin = s * BatchSize;
    IdxType d_begin = d * BatchSize;
    for (long j = 0; j < BatchSize; j++) {
      if (CurrActiveArray[s_begin + j]) {
        uintE newValue = Levels[s_begin + j] + 1;
        if (writeMin(&Levels[d_begin + j], newValue) && CAS(&NextActiveArray[d_begin + j], false, true)) {
          ret = true;
        }
      }
    }
    return ret;
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uintE d) { return cond_true(d); } 
};

template <class vertex>
uintE* Compute_Eval(graph<vertex>& G, std::vector<long> vecQueries, commandLine P) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  uintE* Hops = pbbs::new_array<uintE>(totalNumVertices);
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
    Hops[i] = (uintE)MAXLEVEL;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    Hops[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  parallel_for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);
  while(!Frontier.isEmpty()){
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, BFSLV_F(Hops, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
  return Hops;
}

template <class vertex>
pair<size_t, size_t> Compute_Base(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  intE numVersions = G.batchNum + 1;
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)numVersions;
  fprintf(stderr, "totalNumVertices: %lld\n", totalNumVertices);
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for (long i = 0; i < batch_size; i++) 
  { 
    for (int v = 0; v < numVersions; v++) 
    {
        int idx = vecQueries[i] * (numVersions * batch_size) + v + numVersions * i;
        ShortestPathLen[idx] = 0;
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

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;
  vector<pair<size_t, double>> affinity_tracking;
  vector<pair<size_t, double>> affinity_tracking_one;
  vector<pair<size_t, double>> affinity_tracking_only;
  size_t* overlaps = pbbs::new_array<size_t>(batch_size);

  for (int i = 0; i < batch_size; i++) {
    overlaps[i] = 0;
  }

  vector<double> overlap_scores;
  size_t accumulated_overlap = 0;
  size_t accumulated_overlap_one = 0;
  size_t accumulated_overlap_only= 0;
  size_t peak_activation = 0;
  int peak_iter = 0;
  size_t total_edges = 0;
  // vector<long> frontier_iterations;
  // vector<long> overlapped_iterations;
  // vector<long> accumulated_overlapped_iterations;
  // vector<long> total_activated_iterations;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();
    // cout << "iteration: " << Frontier.size() << endl;
    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size, numVersions), -1, no_dense|remove_duplicates);
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
  // for(int snapshot = 0; snapshot < G.batchNum; snapshot++)
  // {
  //   intE* ShortestPathLen_snapshot = pbbs::new_array<intE>(n);
  //   bool* CurrActiveArray_snapshot = pbbs::new_array<bool>(n);
  //   bool* NextActiveArray_snapshot = pbbs::new_array<bool>(n);
  //   bool* frontier_snapshot = pbbs::new_array<bool>(n);
  //   parallel_for(size_t i = 0; i < n; i++) {
  //     frontier_snapshot[i] = false;
  //   }
  //   for(long i = 0; i < batch_size; i++) {
  //     frontier_snapshot[vecQueries[i]] = true;
  //   }
  //   parallel_for(IdxType i = 0; i < n; i++) {
  //     ShortestPathLen_snapshot[i] = (intE)MAXPATH;
  //     CurrActiveArray_snapshot[i] = false;
  //     NextActiveArray_snapshot[i] = false;
  //   }
  //   ShortestPathLen_snapshot[vecQueries[0]] = 0;
  //   CurrActiveArray_snapshot[vecQueries[0]] = true;

  //   vertexSubset Frontier_snapshot(n, frontier_snapshot);
  //   while(!Frontier_snapshot.isEmpty()){
  //     // mode: no_dense, remove_duplicates (for batch size > 1)
  //     vertexSubset output = edgeMap(G, Frontier_snapshot, DJ_SINGLE_F(ShortestPathLen_snapshot, snapshot, numVersions), -1, no_dense|remove_duplicates);
  //     Frontier_snapshot.del();
  //     Frontier_snapshot = output;

  //     Frontier_snapshot.toDense();
  //     bool* new_d = Frontier_snapshot.d;
  //     Frontier_snapshot.d = nullptr;
  //     vertexSubset Frontier_new(n, new_d);
  //     Frontier_snapshot.del();
  //     Frontier_snapshot = Frontier_new;
  //   }
  //   Frontier_snapshot.del();
  //   intE same = 0;
  //   for(intE node = 0; node < n; node++)
  //   {
  //     if (ShortestPathLen[node * numVersions * batch_size + snapshot * batch_size] == ShortestPathLen_snapshot[node]) {
  //       same++;
  //   }
  //   else
  //   {
  //     cout << "node: " << node << " snapshot: " << snapshot << " " << ShortestPathLen[node * numVersions * batch_size + snapshot * batch_size] << " " << ShortestPathLen_snapshot[node] << endl;
  //   }
  //   }
  //   cout << "snapshot " << snapshot << " same: " << same << endl;
  //   pbbs::delete_array(ShortestPathLen_snapshot, n);
  //   pbbs::delete_array(CurrActiveArray_snapshot, n);
  //   pbbs::delete_array(NextActiveArray_snapshot, n);
  // }

  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(overlaps, batch_size);
  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Separate(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();

  cout << "initializing\n";
  intE** vals = pbbs::new_array<intE*>(batch_size);
  bool** frontiers = pbbs::new_array<bool*>(n);
  for (int i = 0; i < batch_size; i++) {
    vals[i] = pbbs::new_array<intE>(n);
    frontiers[i] = pbbs::new_array<bool>(n);
    parallel_for(IdxType j = 0; j < n; j++) {
      vals[i][j] = (intE)MAXPATH;
      frontiers[i][j] = false;
    }
    vals[i][(IdxType)vecQueries[i]] = 0;
    frontiers[i][(IdxType)vecQueries[i]] = true;
  }

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  bool isConverged = true;
  vertexSubset** vecFs = new vertexSubset*[batch_size];
  // vector<vertexSubset*> vecFs(batch_size);
  for (int i = 0; i < batch_size; i++) {
    // vertexSubset Frontier_new(n, frontiers[i]);

    vecFs[i] = new vertexSubset(n, frontiers[i]);
    // vecFs.push_back(&Frontier_new);
    isConverged = isConverged && vecFs[i]->isEmpty();
  }
  cout << "finished initializing\n";

  while (!isConverged) {
    iteration++;
    // cout << "iteration: " << iteration << endl;
    parallel_for(int j = 0; j < batch_size; j++) {
      vertexSubset output = edgeMap(G, *vecFs[j], DJ_SINGLE_F(vals[j], 0 ,0), -1, no_dense|remove_duplicates);
      vecFs[j]->del();
      *vecFs[j] = output;
    }
    bool isEmpty = true;
    for (int i = 0; i < batch_size; i++) {
      isEmpty = isEmpty && vecFs[i]->isEmpty();
    }
    isConverged = isEmpty;
  }
  // cout << "Total iterations: " << iteration << endl;

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_separate_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, vals[i][j]);
    fclose(fp);
  }
#endif

  // Frontier.del();
  for (int i = 0; i < batch_size; i++) {
    vecFs[i]->del();
    pbbs::delete_array(vals[i], n);
  }

  return make_pair(totalActivated, totalNoOverlap);
}

template <class vertex>
pair<size_t, size_t> Compute_Base_Dynamic(graph<vertex>& G, std::vector<long> vecQueries, queue<long>& queryQueue, commandLine P, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  for(size_t i = 0; i < batch_size; i++) {
    CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
  }

  vertexSubset Frontier(n, frontier);

  
  long iteration = 0;
  size_t totalActivated = 0;
  size_t totalNoOverlap = 0;

  // for async batching
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
    // cout << "iteration: " << iteration << endl;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

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
    if (queued_index != queuedQueries.size()) { // there are more queries to process
      vector<long> empty_slots_vec;
      for (int i = 0; i < batch_size; i++) {
        bool is_finished = true;
        long active_cnt = 0;
        parallel_for(size_t index = 0; index < n; index++) {
          if (CurrActiveArray[index * batch_size + i]) {
            // is_finished = false;
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
            #ifdef OUTPUT 
              long start = vecQueries[index_in_batch];
              char outFileName[300];
              sprintf(outFileName, "SSSP_base_async_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
              FILE *fp;
              fp = fopen(outFileName, "w");
              for (long j = 0; j < n; j++)
                fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + index_in_batch]);
              fclose(fp);
            #endif
            vecQueries[index_in_batch] = nextQ;
            finish_status[index_in_batch] = true;

            CurrActiveArray[nextQ * batch_size + index_in_batch] = true;
            parallel_for(size_t index = 0; index < n; index++) {
              ShortestPathLen[index * batch_size + index_in_batch] = (intE)MAXPATH;
            }
            ShortestPathLen[nextQ * batch_size + index_in_batch] = 0;

            Frontier.toDense();
            bool* new_d = pbbs::new_array<bool>(n);
            parallel_for(size_t ii = 0; ii < n; ii++) {
              new_d[ii] = Frontier.d[ii];
            }
            new_d[nextQ] = true;
            vertexSubset Frontier_new(n, new_d);
            Frontier.del();
            Frontier = Frontier_new;

            // cout << "iteration: " << iteration << ": inserting " << nextQ << " into the batch\n";
            // cout << "current empty slots: " << empty_slots << endl;
          }
        }
      }
    }
    // if (queued_index >= queuedQueries.size() && lastQueries.size() == 0) {
    //   for (int i = 0; i < batch_size; i++) {

    //   }
    // }
  }

  // profiling
  if (should_profile) {
    
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_base_async_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  return make_pair(totalActivated, totalNoOverlap);
}

// This function is modified to compute for batch of queries with different versions
// Value Array now is [node0_{version 0}, node0_{version 1}, node0_{version 2}, node1_{version 0}, node1_{version 1}, node1_{version 2}, ...]
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

    std::cout << "\n=== Processing version subset: ";
    for (intE i = 0; i < versionNum; i++) std::cout << versionNumArray[i] << " ";
    std::cout << "===\n";
#endif

    size_t n = G.n;
    size_t edge_count = G.m;
    long batch_size = vecQueries.size();
#ifndef VERSIONED
    IdxType totalNumVertices = (IdxType)n * (IdxType)versionNum;
#else
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;
#endif
    intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      ShortestPathLen[i] = (intE)MAXPATH;
    }

#ifndef VERSIONED
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * versionNum + v);
        ShortestPathLen[idx] = 0;
      }
    }
#else
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        ShortestPathLen[idx] = 0;
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
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, batch_size), -1, no_dense | remove_duplicates);
#else
      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(ShortestPathLen, batch_size, totalVersion, versionNum, versionNumArray),
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
    std::cout << "Propagation time for version subset: " << t_propagate.totalTime << " seconds" << std::endl;
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
          tmpVal[i] = (intE)MAXPATH;
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;
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
        // intE error = 0;
        // for (intE node = 0; node < n; node++) {
        //   IdxType idx = ((IdxType)node * batch_size + query) * versionNum + v;
        //   if (ShortestPathLen[idx] != tmpVal[node]) {
        //     error++;
        //   }
        // }
        // std::cout << "query: " << query << " version: " << version << " error: " << error << std::endl;
        pbbs::delete_array(tmpVal, n);
      }
    }
    
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(ShortestPathLen, totalNumVertices);
  }
#ifdef VERSIONED
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
  std::cout << "Total baseline time: " << baseTime << " seconds\n";
  std::cout << "Speedup: " << baseTime / total_propagation_time << std::endl;
#endif
  return make_pair(0, 0); // totalActivated, totalNoOverlap
}

template <class vertex>
pair<size_t, size_t> Compute_Delay(graph<vertex>& G, std::vector<long> vecQueries, commandLine P, std::vector<int> defer_vec, bool should_profile) {
  size_t n = G.n;
  size_t edge_count = G.m;
  long batch_size = vecQueries.size();
  IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size;
  intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
  bool* CurrActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* NextActiveArray = pbbs::new_array<bool>(totalNumVertices);
  bool* frontier = pbbs::new_array<bool>(n);
  parallel_for(size_t i = 0; i < n; i++) {
    frontier[i] = false;
  }
  // for delaying initialization
  for(long i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      frontier[vecQueries[i]] = true;
    }
  }
  // for(long i = 0; i < batch_size; i++) {
  //   frontier[vecQueries[i]] = true;
  // }
  parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
    ShortestPathLen[i] = (intE)MAXPATH;
    CurrActiveArray[i] = false;
    NextActiveArray[i] = false;
  }
  for(long i = 0; i < batch_size; i++) {
    ShortestPathLen[(IdxType)batch_size * (IdxType)vecQueries[i] + (IdxType)i] = 0;
  }
  for(size_t i = 0; i < batch_size; i++) {
    if (defer_vec[i] == 0) {
      CurrActiveArray[(IdxType)vecQueries[i] * (IdxType)batch_size + (IdxType)i] = true;
    }
  }

  vertexSubset Frontier(n, frontier);

  // for profiling
  long iteration = 0;
  size_t totalActivated = 0;

  while(!Frontier.isEmpty()){
    iteration++;
    totalActivated += Frontier.size();

    // mode: no_dense, remove_duplicates (for batch size > 1)
    vertexSubset output = edgeMap(G, Frontier, DJ_F(ShortestPathLen, CurrActiveArray, NextActiveArray, batch_size), -1, no_dense|remove_duplicates);

    Frontier.del();
    Frontier = output;

    // Checking for delayed queries.
    // TODO: possibly optimization.
    // timer t_checking;
    // t_checking.start();
    // Frontier.toDense();
    // bool* new_d = pbbs::new_array<bool>(n);
    // parallel_for(size_t i = 0; i < n; i++) {
    //   new_d[i] = Frontier.d[i];
    // }

    
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
    // t_checking.stop();
    // t_checking.reportTotal("checking delaying");

    std::swap(CurrActiveArray, NextActiveArray);
    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      NextActiveArray[i] = false;
    }
  }

#ifdef OUTPUT 
  for (int i = 0; i < batch_size; i++) {
    long start = vecQueries[i];
    char outFileName[300];
    sprintf(outFileName, "SSSP_delay_output_src%ld.%ld.%ld.out", start, edge_count, batch_size);
    FILE *fp;
    fp = fopen(outFileName, "w");
    for (long j = 0; j < n; j++)
      fprintf(fp, "%ld %d\n", j, ShortestPathLen[j * batch_size + i]);
    fclose(fp);
  }
#endif

  Frontier.del();
  pbbs::delete_array(CurrActiveArray, totalNumVertices);
  pbbs::delete_array(NextActiveArray, totalNumVertices);
  pbbs::delete_array(ShortestPathLen, totalNumVertices);
  return make_pair(totalActivated, 0);
}

template <class vertex>
pair<size_t, size_t> Compute_Delay_Skipping(graph<vertex>& G,
  std::vector<long> vecQueries, commandLine P,
  std::vector<int> defer_vec, bool should_profile)
{
#ifdef VERSIONED
  // auto totalVersion = G.batchNum + 1;
  auto totalVersion = G.batchNum;
  auto parallelVersionNum = P.getOptionIntValue("-parallelVersion", 4);
  double total_propagation_time = 0.0;

  for (intE versionStart = 0; versionStart < totalVersion; versionStart += parallelVersionNum) {
    intE versionNum = std::min<intE>(parallelVersionNum, totalVersion - versionStart);
    intE* versionNumArray = newA(intE, versionNum);
    for (intE i = 0; i < versionNum; i++) {
      versionNumArray[i] = versionStart + i;
    }

    std::cout << "\n=== Processing version subset: ";
    for (intE i = 0; i < versionNum; i++) std::cout << versionNumArray[i] << " ";
    std::cout << "===\n";
#endif

    size_t n = G.n;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
      ShortestPathLen[i] = (intE)MAXPATH;
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        ShortestPathLen[idx] = 0;
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
        DJ_SKIP_F(ShortestPathLen, batch_size, totalVersion, versionNum, versionNumArray),
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
    std::cout << "Propagation time for version subset: " << t_propagate.totalTime << " seconds" << std::endl;
    total_propagation_time += t_propagate.totalTime;
#endif

    Frontier.del();
#ifdef VERSIONED
    pbbs::delete_array(versionNumArray, versionNum);
#endif
    pbbs::delete_array(ShortestPathLen, totalNumVertices);
  }
#ifdef VERSIONED
  std::cout << "\nTotal propagation time (excluding baseline): " << total_propagation_time << " seconds\n";
#endif
  return make_pair(0, 0);
}

// This function is modified to compute for batch of queries with different versions
// Value Array now is [node0_{version 0}, node0_{version 1}, node0_{version 2}, node1_{version 0}, node1_{version 1}, node1_{version 2}, ...]
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

    std::cout << "\n=== Processing version subset: ";
    for (intE i = 0; i < versionNum; i++) std::cout << versionNumArray[i] << " ";
    std::cout << "===\n";
#endif

    size_t n = G.n;
    size_t edge_count = G.m;
    long batch_size = vecQueries.size();
#ifndef VERSIONED
    IdxType totalNumVertices = (IdxType)n * (IdxType)versionNum;
#else
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;
#endif
    intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
    bool* frontier = pbbs::new_array<bool>(n);

    parallel_for(size_t i = 0; i < n; i++) frontier[i] = false;
    for (long i = 0; i < batch_size; i++) frontier[vecQueries[i]] = true;

    parallel_for(IdxType i = 0; i < totalNumVertices; i++) {
      ShortestPathLen[i] = (intE)MAXPATH;
    }

#ifndef VERSIONED
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * versionNum + v);
        ShortestPathLen[idx] = 0;
      }
    }
#else
    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        ShortestPathLen[idx] = 0;
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
      vertexSubset output = edgeMap(G, Frontier, DJ_SKIP_F(ShortestPathLen, batch_size), -1, no_dense | remove_duplicates);
#else
      vertexSubset output = edgeMap(G, Frontier,
        DJ_SKIP_F(ShortestPathLen, batch_size, totalVersion, versionNum, versionNumArray),
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
    std::cout << "Propagation time for version subset: " << t_propagate.totalTime << " seconds" << std::endl;
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
          tmpVal[i] = (intE)MAXPATH;
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;
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
    pbbs::delete_array(ShortestPathLen, totalNumVertices);
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
  // auto totalVersion = G.batchNum + 1;
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

    std::cout << "\n=== Processing version subset: ";
    for (intE i = 0; i < versionNum; i++) std::cout << versionNumArray[i] << " ";
    std::cout << "===\n";
#endif

    size_t n = G.n;
    long batch_size = vecQueries.size();
    IdxType totalNumVertices = (IdxType)n * (IdxType)batch_size * (IdxType)versionNum;

    intE* ShortestPathLen = pbbs::new_array<intE>(totalNumVertices);
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
      ShortestPathLen[i] = (intE)MAXPATH;
    }

    for (long i = 0; i < batch_size; i++) {
      uintE src = vecQueries[i];
      for (intE v = 0; v < versionNum; v++) {
        IdxType idx = ((IdxType)src * batch_size + i) * versionNum + v;
        ShortestPathLen[idx] = 0;
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
        DJ_SKIP_F(ShortestPathLen, batch_size, totalVersion, versionNum, versionNumArray),
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
    std::cout << "Propagation time for version subset: " << t_propagate.totalTime << " seconds" << std::endl;
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
          tmpVal[i] = (intE)MAXPATH;
          tmpFrontier[i] = false;
        }

        tmpVal[src] = 0;
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
    pbbs::delete_array(ShortestPathLen, totalNumVertices);
  }
#ifdef VERSIONED
  return make_pair(total_propagation_time, baseTime);
#else
  return make_pair(0.0, 0.0);
#endif  
}
