#ifndef GRID_UTILS_H
#define GRID_UTILS_H

inline int 
get_nb_count(char**& sts, int x, int y)
{
  int cnt = 0;
  if (sts[y][x-1]>0)
    cnt++;
  if (sts[y-1][x]>0)
    cnt++;
  if (sts[y][x+1]>0)
    cnt++;
  if (sts[y+1][x]>0)
    cnt++;
  return cnt;
}

inline int
get_nb_corners_count(char**& sts, int x, int y)
{
  int cnt = 0;
  if (sts[y-1][x-1]>0)
    cnt++;
  if (sts[y-1][x+1]>0)
    cnt++;
  if (sts[y+1][x-1]>0)
    cnt++;
  if (sts[y+1][x+1]>0)
    cnt++;
  return cnt;
}

inline int
get_nb8_count(char**& sts, int x, int y)
{
  int cnt = 0;
  cnt = get_nb_count(sts, x, y);
  cnt += get_nb_corners_count(sts, x, y);

  return cnt;
}

inline int
get_corner_nb_count(char**& status, int x, int y)
{
  int cnt = 0;
  if (status[y-1][x-1]>0)
    cnt++;
  if (status[y-1][x+1]>0)
    cnt++;
  if (status[y+1][x-1]>0)
    cnt++;
  if (status[y+1][x+1]>0)
    cnt++;
  return cnt;
}

inline int
diff_neighbors_8(int x, int y)
{
  int diff4[5] = {4, 2, 0, -2, -4};
  int nb4 = get_nb_count(status, x, y);

  int r = diff4[nb4];

  int nb_corners = get_corner_nb_count(status, x,y);
  // seems to work well: 
  int diff_corners[5] = {4, 2, 0, -2, -4};  
  //int diff_corners[5] = {2, 1.5, 1, 0.5, 0};  
  r += diff_corners[nb_corners];
  return r;
}

#endif //GRID_UTILS_H
