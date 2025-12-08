//std::numeric_limits<float>::infinity();
#include <iostream>
#include <queue>
#include <vector>
#include <functional>
#include <math.h>
#include <chrono>
#include <random>
#include <tuple>
#include <limits>
#include <algorithm>
#include <fstream>
#include <omp.h>
#include <immintrin.h>
#include <atomic>

template<int D = 4>
class DaryHeap {
private:
    std::vector<int> indices;
    std::vector<float> f_costs;
    std::vector<float> g_costs;
    std::vector<int> parents;

    void siftUp(int pos) {
        int parent_pos = (pos - 1) / D;
        float my_f = f_costs[pos];
        int my_idx = indices[pos];
        float my_g = g_costs[pos];
        int my_parent = parents[pos];

        while(pos > 0 && my_f < f_costs[parent_pos]) {
            indices[pos] = indices[parent_pos];
            f_costs[pos] = f_costs[parent_pos];
            g_costs[pos] = g_costs[parent_pos];
            parents[pos] = parents[parent_pos];
            pos = parent_pos;
            parent_pos = (pos - 1) / D;
        }

        indices[pos] = my_idx;
        f_costs[pos] = my_f;
        g_costs[pos] = my_g;
        parents[pos] = my_parent;
    }

    void siftDown(int pos) {
        int size = f_costs.size();
        float my_f = f_costs[pos];
        int my_idx = indices[pos];
        float my_g = g_costs[pos];
        int my_parent = parents[pos];

        while(true) {
            int first_child = D * pos + 1;
            if(first_child >= size) break;

            int smallest_child = first_child;
            float smallest_f = f_costs[first_child];

            int last_child = std::min(first_child + D, size);
            for(int child = first_child + 1; child < last_child; child++) {
                if(f_costs[child] < smallest_f) {
                    smallest_child = child;
                    smallest_f = f_costs[child];
                }
            }

            if(my_f <= smallest_f) break;

            indices[pos] = indices[smallest_child];
            f_costs[pos] = f_costs[smallest_child];
            g_costs[pos] = g_costs[smallest_child];
            parents[pos] = parents[smallest_child];
            pos = smallest_child;
        }

        indices[pos] = my_idx;
        f_costs[pos] = my_f;
        g_costs[pos] = my_g;
        parents[pos] = my_parent;
    }

public:
    void reserve(size_t capacity) {
        indices.reserve(capacity);
        f_costs.reserve(capacity);
        g_costs.reserve(capacity);
        parents.reserve(capacity);
    }

    void push(int idx, float f, float g, int parent) {
        indices.push_back(idx);
        f_costs.push_back(f);
        g_costs.push_back(g);
        parents.push_back(parent);
        siftUp(indices.size() - 1);
    }

    std::tuple<int, float, float, int> pop() {
        int result_idx = indices[0];
        float result_f = f_costs[0];
        float result_g = g_costs[0];
        int result_parent = parents[0];

        int last = indices.size() - 1;
        indices[0] = indices[last];
        f_costs[0] = f_costs[last];
        g_costs[0] = g_costs[last];
        parents[0] = parents[last];

        indices.pop_back();
        f_costs.pop_back();
        g_costs.pop_back();
        parents.pop_back();

        if(!indices.empty())
            siftDown(0);

        return {result_idx, result_f, result_g, result_parent};
    }

    bool empty() const { return indices.empty(); }
    size_t size() const { return indices.size(); }
};

class Graf {
private:
  const std::vector<float>&cvorovi_v2;
  int broj_redova, broj_kolona;

  static float Heuristika(int current_x, int current_y, int goal_x, int goal_y, float w_min) {
    int dx = std::abs(goal_x - current_x);
    int dy = std::abs(goal_y - current_y);
    return w_min * (float)std::max(dx, dy);
  }

public:
  Graf(const std::vector<float> &c, int r, int k): cvorovi_v2(c) {
    broj_redova = r;
    broj_kolona = k;
  }

  // Osnovna verzija - za poređenje
  std::pair<float, std::vector<int>> A_star_v3(std::pair<int,int>start, std::pair<int,int>finish) {
    if(start.first < 0 || start.first >= broj_redova ||
       start.second < 0 || start.second >= broj_kolona ||
       finish.first < 0 || finish.first >= broj_redova ||
       finish.second < 0 || finish.second >= broj_kolona) {
        throw std::domain_error("Start ili finis van granica");
    }

    int index_trenutniCvor = start.first * broj_kolona + start.second;
    int index_finishCvor = finish.first * broj_kolona + finish.second;

    const int dx[8] = { 0, -1, 0, 1, -1, -1, 1, 1 };
    const int dy[8] = { 1, 0, -1, 0, 1, -1, 1, -1 };

    std::vector<std::tuple<bool, float, int>> obradjeniCvorovi(broj_redova * broj_kolona, std::make_tuple(false, std::numeric_limits<float>::infinity(), 0));
    DaryHeap<4> minHeap;
    minHeap.reserve(broj_redova * broj_kolona / 20);
    minHeap.push(index_trenutniCvor, 0, 0, 0);

    std::vector<float> heuristika_za_okolne_cvorove(8);
    std::tuple<int,float,float,int> trenutniCvor;

    while(!minHeap.empty()) {
        trenutniCvor = minHeap.pop();

        if(std::get<0>(trenutniCvor) == index_finishCvor){
          std::vector<int> putanja;
          index_trenutniCvor = std::get<0>(trenutniCvor);
          obradjeniCvorovi[index_trenutniCvor] = {true, std::get<2>(trenutniCvor), std::get<3>(trenutniCvor)};

          int trenutni = index_finishCvor;
          int start_index = start.first * broj_kolona + start.second;

          while(trenutni != start_index) {
              putanja.push_back(trenutni);
              trenutni = std::get<2>(obradjeniCvorovi[trenutni]);
          }
          putanja.push_back(start_index);

          return {std::get<2>(trenutniCvor), putanja};
        }

        index_trenutniCvor = std::get<0>(trenutniCvor);
        if(std::get<0>(obradjeniCvorovi[index_trenutniCvor]))
            continue;

        int x = index_trenutniCvor / broj_kolona;
        int y = index_trenutniCvor % broj_kolona;

        for(int i = 0; i < 8; i++)
            heuristika_za_okolne_cvorove[i] = Heuristika(x + dx[i], y + dy[i], finish.first, finish.second, 1.0f);

        for(int i = 0; i < 8; i++) {
            int nx = x + dx[i];
            int ny = y + dy[i];
            if(nx < 0 || nx >= broj_redova || ny < 0 || ny >= broj_kolona)
                continue;

            int weight_index = index_trenutniCvor * 8 + i;
            if(cvorovi_v2[weight_index] == std::numeric_limits<float>::infinity())
                continue;

            int index_susjeda = nx * broj_kolona + ny;
            if(!std::get<0>(obradjeniCvorovi[index_susjeda])) {
                float g_novo = std::get<2>(trenutniCvor) + cvorovi_v2[weight_index];
                float f_novo = g_novo + heuristika_za_okolne_cvorove[i];
                if(f_novo < std::get<1>(obradjeniCvorovi[index_susjeda])){
                  minHeap.push(index_susjeda, f_novo, g_novo, index_trenutniCvor);
                  std::get<1>(obradjeniCvorovi[index_susjeda]) = f_novo;
                }
            }
        }

        obradjeniCvorovi[index_trenutniCvor] = {true, std::get<2>(trenutniCvor), std::get<3>(trenutniCvor)};
    }

    return {-1., std::vector<int>(0)};
  }

  std::pair<float, std::vector<int>> A_star_hierarchical(std::pair<int,int>start, std::pair<int,int>finish) {
    if(start.first < 0 || start.first >= broj_redova ||
       start.second < 0 || start.second >= broj_kolona ||
       finish.first < 0 || finish.first >= broj_redova ||
       finish.second < 0 || finish.second >= broj_kolona) {
        throw std::domain_error("Start ili finis van granica");
    }

    // Podjela na 4 regiona (2x2)
    int mid_row = broj_redova / 2;
    int mid_col = broj_kolona / 2;

    std::cout << "Grid podijeljen na 4 regiona sa granicom: row=" << mid_row << ", col=" << mid_col << std::endl;

    // KORAK 1: Nađi presjecanja prave linije (start->finish) sa granicama regiona
    std::vector<std::pair<int,int>> waypoints;
    waypoints.push_back(start);  // Dodaj start

    // Parametarska jednačina prave: P(t) = start + t*(finish - start), t ∈ [0,1]
    int dx = finish.first - start.first;
    int dy = finish.second - start.second;

    std::cout << "Prava linija: od (" << start.first << "," << start.second
              << ") do (" << finish.first << "," << finish.second << ")" << std::endl;

    // Provjeri presjecanje sa horizontalnom granicom (mid_row)
    if((start.first < mid_row && finish.first >= mid_row) ||
       (start.first >= mid_row && finish.first < mid_row)) {
        // Prava linija presijecа horizontalnu granicu
        // Nađi t gdje je x = mid_row: start.first + t*dx = mid_row
        float t = (float)(mid_row - start.first) / (float)dx;
        int intersect_col = start.second + (int)(t * dy);

        // Provjeri da li je presjecanje unutar grida
        if(intersect_col >= 0 && intersect_col < broj_kolona) {
            // Dodaj obje strane granice
            if(start.first < mid_row) {
                waypoints.push_back({mid_row - 1, intersect_col});  // Zadnji red gornjeg
                waypoints.push_back({mid_row, intersect_col});      // Prvi red donjeg
            } else {
                waypoints.push_back({mid_row, intersect_col});      // Prvi red donjeg
                waypoints.push_back({mid_row - 1, intersect_col});  // Zadnji red gornjeg
            }
            std::cout << "Horizontalno presjecanje na: (" << mid_row-1 << "," << intersect_col
                      << ") i (" << mid_row << "," << intersect_col << ")" << std::endl;
        }
    }

    // Provjeri presjecanje sa vertikalnom granicom (mid_col)
    if((start.second < mid_col && finish.second >= mid_col) ||
       (start.second >= mid_col && finish.second < mid_col)) {
        // Prava linija presijecа vertikalnu granicu
        // Nađi t gdje je y = mid_col: start.second + t*dy = mid_col
        float t = (float)(mid_col - start.second) / (float)dy;
        int intersect_row = start.first + (int)(t * dx);

        // Provjeri da li je presjecanje unutar grida
        if(intersect_row >= 0 && intersect_row < broj_redova) {
            // Dodaj obje strane granice
            if(start.second < mid_col) {
                waypoints.push_back({intersect_row, mid_col - 1}); // Zadnja kolona lijeve
                waypoints.push_back({intersect_row, mid_col});     // Prva kolona desne
            } else {
                waypoints.push_back({intersect_row, mid_col});     // Prva kolona desne
                waypoints.push_back({intersect_row, mid_col - 1}); // Zadnja kolona lijeve
            }
            std::cout << "Vertikalno presjecanje na: (" << intersect_row << "," << mid_col-1
                      << ") i (" << intersect_row << "," << mid_col << ")" << std::endl;
        }
    }

    waypoints.push_back(finish);  // Dodaj finish

    std::cout << "Ukupno waypoints: " << waypoints.size() << std::endl;
    for(int i = 0; i < waypoints.size(); i++) {
        std::cout << "  Waypoint " << i << ": (" << waypoints[i].first << "," << waypoints[i].second << ")" << std::endl;
    }

    // KORAK 2: Paralelno traži puteve između uzastopnih waypoints
    int num_segments = waypoints.size() - 1;
    std::vector<std::pair<float, std::vector<int>>> segment_paths(num_segments);

    std::cout << "\nTrazi " << num_segments << " segmenata paralelno..." << std::endl;

    #pragma omp parallel for schedule(dynamic) num_threads(4)
    for(int i = 0; i < num_segments; i++) {
        auto seg_start = waypoints[i];
        auto seg_end = waypoints[i + 1];

        #pragma omp critical
        {
            std::cout << "Thread " << omp_get_thread_num() << " trazi segment " << i
                      << ": (" << seg_start.first << "," << seg_start.second << ") -> ("
                      << seg_end.first << "," << seg_end.second << ")" << std::endl;
        }

        segment_paths[i] = A_star_v3(seg_start, seg_end);

        #pragma omp critical
        {
            if(segment_paths[i].first >= 0) {
                std::cout << "Thread " << omp_get_thread_num() << " zavrsen segment " << i
                          << " sa cijenom " << segment_paths[i].first << std::endl;
            } else {
                std::cout << "Thread " << omp_get_thread_num() << " NIJE NASAO PUT za segment " << i << std::endl;
            }
        }
    }

    // KORAK 3: Provjeri da li su svi segmenti pronađeni
    float total_cost = 0.0f;
    for(int i = 0; i < num_segments; i++) {
        if(segment_paths[i].first < 0) {
            std::cout << "GRESKA: Segment " << i << " nije pronadjen!" << std::endl;
            return {-1., std::vector<int>(0)};
        }
        total_cost += segment_paths[i].first;
    }

    // KORAK 4: Spoji sve segmente u jedan put
    std::vector<int> full_path;

    for(int i = 0; i < num_segments; i++) {
        const auto& seg_path = segment_paths[i].second;

        if(i == 0) {
            // Prvi segment - dodaj sve
            full_path.insert(full_path.end(), seg_path.begin(), seg_path.end());
        } else {
            // Ostali segmenti - preskoči prvu tačku (duplikat sa prethodnim)
            full_path.insert(full_path.end(), seg_path.begin() + 1, seg_path.end());
        }
    }

    std::cout << "\nFinalni put spojen! Ukupna cijena: " << total_cost << std::endl;
    std::cout << "Ukupno cvorova u putu: " << full_path.size() << std::endl;

    return {total_cost, full_path};
  }
};

int main()
{
    const int rows = 10000;
    const int cols = 10000;
    std::pair<int,int> start = {0, 0};
    std::pair<int,int> finish = {rows - 1, cols - 1};

    std::random_device rd;
    std::mt19937 gen(42);  // Fixed seed
    std::uniform_real_distribution<float> dis(1.0f, 10.0f);
    std::uniform_real_distribution<float> wall(0.0f, 1.0f);

    std::cout << "Priprema mreze " << rows << "x" << cols << "..." << std::endl;

    std::vector<float> cvorovi_flat(rows * cols * 8, std::numeric_limits<float>::infinity());

    const int dx[8]  = { 0, -1, 0, 1, -1, -1, 1, 1 };
    const int dy[8]  = { 1,  0,-1, 0,  1, -1, 1,-1 };
    const int opp[8] = { 2,  3, 0, 1,  7,  6, 5, 4 };

    for (int x = 0; x < rows; x++) {
      for (int y = 0; y < cols; y++) {
        int u = x * cols + y;
        for (int i = 0; i < 8; i++) {
          int nx = x + dx[i];
          int ny = y + dy[i];
          if (nx < 0 || nx >= rows || ny < 0 || ny >= cols) continue;
          int v = nx * cols + ny;
          if (u < v) {
            float w;
            if (wall(gen) < 0.05f) w = std::numeric_limits<float>::infinity();
            else w = dis(gen);
            cvorovi_flat[u * 8 + i] = w;
            cvorovi_flat[v * 8 + opp[i]] = w;
          }
        }
      }
    }

    Graf graf(cvorovi_flat, rows, cols);

    std::cout << "\n========================================" << std::endl;
    std::cout << "=== STANDARDNI A* (v3) ===" << std::endl;
    std::cout << "========================================" << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto [cost_v3, path_v3] = graf.A_star_v3(start, finish);
    auto t2 = std::chrono::high_resolution_clock::now();

    if(cost_v3 >= 0) {
        std::cout << "Put pronadjen!" << std::endl;
        std::cout << "Cijena: " << cost_v3 << std::endl;
        std::cout << "Duzina puta: " << path_v3.size() << " cvorova" << std::endl;
    } else {
        std::cout << "Put nije pronaden!" << std::endl;
    }

    double time_v3 = std::chrono::duration<double>(t2 - t1).count();
    std::cout << "Vrijeme izvrsavanja: " << time_v3 << " sekundi" << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "=== HIJERARHIJSKI A* (4 REGIONA) ===" << std::endl;
    std::cout << "========================================" << std::endl;

    auto t3 = std::chrono::high_resolution_clock::now();
    auto [cost_hier, path_hier] = graf.A_star_hierarchical(start, finish);
    auto t4 = std::chrono::high_resolution_clock::now();

    if(cost_hier >= 0) {
        std::cout << "Put pronadjen!" << std::endl;
        std::cout << "Cijena: " << cost_hier << std::endl;
        std::cout << "Duzina puta: " << path_hier.size() << " cvorova" << std::endl;
    } else {
        std::cout << "Put nije pronaden!" << std::endl;
    }

    double time_hier = std::chrono::duration<double>(t4 - t3).count();
    std::cout << "Vrijeme izvrsavanja: " << time_hier << " sekundi" << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "=== POREDENJE ===" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Standardni A* (v3): " << time_v3 << " s" << std::endl;
    std::cout << "Hijerarhijski A*:   " << time_hier << " s" << std::endl;

    if(time_v3 > time_hier) {
        std::cout << "\nSpeedup: " << (time_v3 / time_hier) << "x brze!" << std::endl;
    } else {
        std::cout << "\nSlowdown: " << (time_hier / time_v3) << "x sporije" << std::endl;
    }

    if(cost_v3 >= 0 && cost_hier >= 0) {
        float diff = std::abs(cost_v3 - cost_hier);
        float percent = (diff / cost_v3) * 100.0f;
        std::cout << "Razlika u cijeni: " << diff << " (" << percent << "%)" << std::endl;
    }

    return 0;
}
