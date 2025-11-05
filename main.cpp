#include <iostream>
#include <queue>
#include <vector>
#include <functional>
#include <math.h>
//std::numeric_limits<float>::infinity();

class Graf {
private:
  const std::vector<std::vector<float>> &cvorovi;
  int broj_redova, broj_kolona;

  static float Heuristika(int current_x, int current_y, int goal_x, int goal_y, float w_min) {
    int dx = std::abs(goal_x - current_x);
    int dy = std::abs(goal_y - current_y);
    return w_min * std::sqrt(dx * dx + dy * dy);
  }
public:
  Graf(const std::vector<std::vector<float>> &c, int r, int k) : cvorovi(c) {
    broj_redova = r;
    broj_kolona = k;
  }
  //Osnovna sekvencijalna metoda
  float A_star(std::pair<int,int>start, std::pair<int,int>finish) {
    int index_trenutniCvor = start.first*broj_kolona+start.second;
    int index_finishCvor = finish.first*broj_kolona+finish.second;
    const int dx[8] = { 1, 0, -1,  0, 1, -1,  1, -1 };
    const int dy[8] = { 0, -1,  0, 1, -1,  -1, 1, 1 };
    std::vector<std::tuple<bool, float, int>>/*obradjen, g(n), od kojeg cvora je dobio */ obradjeniCvorovi(broj_redova*broj_kolona, std::make_tuple(false, 0.0f, 0));

    auto cmp = [](const std::tuple<int,float, float,int> &a, const std::tuple<int,float, float,int> &b){
        return std::get<1>(a) > std::get<1>(b);
    };
    std::priority_queue<std::tuple<int,float,float,int>/*cvor, f(n), g(n),od_koga*/,std::vector<std::tuple<int,float,float,int>>,decltype(cmp)> minHeap(cmp);

    minHeap.push({index_trenutniCvor, 0,0,0});
    std::vector<float> heuristika_za_okolne_cvorove(8);

    while(minHeap.size() > 0) {
      std::tuple<int,float,float,int> trenutniCvor = minHeap.top();
      minHeap.pop();
      if(std::get<0>(trenutniCvor)==index_finishCvor){
        return std::get<2>(trenutniCvor);
      }
      if(std::get<0>(obradjeniCvorovi[std::get<0>(trenutniCvor)]) == true){
         continue;
      }
      index_trenutniCvor = std::get<0>(trenutniCvor);
      int x = index_trenutniCvor / broj_kolona;
      int y = index_trenutniCvor % broj_kolona;

      heuristika_za_okolne_cvorove[0] = Heuristika(x+1,y, finish.first, finish.second,1.);
      heuristika_za_okolne_cvorove[1] = Heuristika(x,y-1, finish.first, finish.second,1.);
      heuristika_za_okolne_cvorove[2] = Heuristika(x-1,y, finish.first, finish.second,1.);
      heuristika_za_okolne_cvorove[3] = Heuristika(x,y+1, finish.first, finish.second,1.);
      heuristika_za_okolne_cvorove[4] = Heuristika(x+1,y-1, finish.first, finish.second,1.);//Gore desno
      heuristika_za_okolne_cvorove[5] = Heuristika(x-1,y-1, finish.first, finish.second,1.);//Gore lijevo
      heuristika_za_okolne_cvorove[6] = Heuristika(x+1,y+1, finish.first, finish.second,1.);//Dole desno
      heuristika_za_okolne_cvorove[7] = Heuristika(x-1,y+1, finish.first, finish.second,1.);//Dole lijevo

      for(int i = 0; i < 8; i++) {
        int nx = x + dx[i];
        int ny = y + dy[i];
        int index_susjeda = nx * broj_kolona + ny;
        minHeap.push({index_susjeda,heuristika_za_okolne_cvorove[i] + std::get<2>(trenutniCvor)+cvorovi[std::get<0>(trenutniCvor)][i], std::get<2>(trenutniCvor)+cvorovi[std::get<0>(trenutniCvor)][i], std::get<0>(trenutniCvor)});
      }
      obradjeniCvorovi[std::get<0>(trenutniCvor)] = {true, std::get<2>(trenutniCvor), std::get<3>(trenutniCvor)};
    }

  }
};
int main()
{
  const int rows = 10;
  const int cols = 10;

  std::vector<std::vector<float>> cvorovi(rows*cols, std::vector<float>(8, 1.0f));

  cvorovi[21][6] = 150.0f;
  //cvorovi[11][4] = 3.0f;

  std::pair<int,int> start = {2,1};
  std::pair<int,int> finish = {4,3};

  Graf graf(cvorovi, rows, cols);
  float path_cost = graf.A_star(start, finish);

  if(path_cost >= 0)
      std::cout << "Najkraci put od (" << start.first << "," << start.second << ") do ("
                << finish.first << "," << finish.second << ") ima cijenu: " << path_cost << "\n";
  return 0;
}
