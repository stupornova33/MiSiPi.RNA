#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
DataFrame mergePileups(std::vector<int> start, std::vector<int> end, std::vector<int> pileups_pos, std::vector<int> pileups_count){
   int startSize = start.size();
   int pileupsSize = pileups_pos.size();

   int memory_reserve = startSize * 2;
   std::vector<int> pos_results;
   pos_results.reserve(memory_reserve);
   std::vector<int> count_results;
   count_results.reserve(memory_reserve);

   //iterate through the start position in all reads
   for (int i = 0; i < startSize; i++) {
      int start_coord = start[i];
      int end_coord = end[i];

      //std::vector<int> vec;

      //iterate through all the positions in one read
      for (int j = start_coord; j <= end_coord; j++) {
         //vec.push_back(j);
         bool found = false;
         for (int k = 0; k < pileupsSize; k++) {
            if (j == pileups_pos[k]) {
               pos_results.emplace_back(j);
               count_results.emplace_back(pileups_count[k]);
               found = true;
            }
         }
         if (found == false) {
            pos_results.emplace_back(j);
            count_results.emplace_back(0);
         }
      }
   }

   DataFrame df = DataFrame::create(Named("pos") = pos_results, Named("count") = count_results);
   return df;
}

//[[Rcpp::export]]
DataFrame get_nearby(IntegerVector f_start, IntegerVector f_end, IntegerVector r_start, IntegerVector r_end, int max_dist, int dfsize){
// Calculate the distance between overlapping read pairs

   int memory_reserve = dfsize;

   std::vector<int> res_start_r1;
   res_start_r1.reserve(memory_reserve);

   std::vector<int> res_end_r1;
   res_end_r1.reserve(memory_reserve);

   std::vector<int> res_wx;
   res_wx.reserve(memory_reserve);

   std::vector<int> res_start_r2;
   res_start_r2.reserve(memory_reserve);

   std::vector<int> res_end_r2;
   res_end_r2.reserve(memory_reserve);

   std::vector<int> res_wy;
   res_wy.reserve(memory_reserve);

   std::vector<int> res_dist;
   res_dist.reserve(memory_reserve);





   for(int i = 0; i < dfsize; i++){
      int r1_start = f_start[i];
      int r1_end = f_end[i];
      int r2_start = r_start[i];
      int r2_end = r_end[i];
      int wx = r1_end - r1_start + 1;
      int wy = r2_end - r2_start + 1;

      // compute the distance between end of r1 and start of r2
      int dist = r2_start - r1_end;

      // matched sequences cannot be overlapping (from same strand)

      if((r1_end < r2_start) && (dist <= max_dist) && (dist > 0))
      {
         res_start_r1.emplace_back(r1_start);
         res_end_r1.emplace_back(r1_end);
         res_wx.emplace_back(wx);
         res_start_r2.emplace_back(r2_start);
         res_end_r2.emplace_back(r2_end);
         res_dist.emplace_back(abs(dist));
         res_wy.emplace_back(wy);
      }


   }

   DataFrame df = DataFrame::create(Named("start_r1") = res_start_r1, Named("end_r1") = res_end_r1, Named("widthx") = res_wx, Named("start_r2") = res_start_r2,
                                    Named("end_r2") = res_end_r2,
                                    Named("widthy") = res_wy, Named("dist") = res_dist);
   return df;


}


//[[Rcpp::export]]
DataFrame getFastas(std::string input_str, IntegerVector start, IntegerVector stop, int size){
   int memory_reserve = size;
   std::vector<std::string> fastas;
   fastas.reserve(memory_reserve);
   std::vector<int> res_start;
   res_start.reserve(memory_reserve);
   std::vector<int>res_end;
   res_end.reserve(memory_reserve);

   for (int i = 0; i < size; i++) {
      int st = start[i];
      int end = stop[i];
      res_start.emplace_back(st);
      res_end.emplace_back(end);
      fastas.push_back(input_str.substr(st, end - st + 1));
   }
   DataFrame df = DataFrame::create( Named("start") = res_start, Named("stop") = res_end, Named("Seq") = fastas);
   return(df);
}


//[[Rcpp::export]]
std::vector<std::string> convertU(std::vector<std::string> input_str, int size){
   std::vector<std::string> converted;
   converted.reserve(size);

   for(int i = 0; i < size; i++){
      std::string current_seq = input_str[i];
      std::string converted_str = "";
      int str_size = current_seq.size();
      for(int j = 0; j < str_size; j++){
         if(current_seq[j] == 'T'){
            converted_str.append("U");
         } else {
            converted_str.push_back(current_seq[j]);
         }
      }
      converted.emplace_back(converted_str);
   }
   return converted;
}


//[[Rcpp::export]]

DataFrame getPileups(std::vector<int> dtpos, std::vector<int> dtcount, std::vector<int> start_r1, std::vector<int> end_r1,
                     std::vector<int> start_r2, std::vector<int> end_r2){
//calculates the average pileup for a read using counts at each nucleotide position
//takes position from read pileups df, count at that position, start pos of candidate read, end pos of candidate read

   int memory_reserve = dtpos.size();

   std::vector<int> res_start_r1;
   res_start_r1.reserve(memory_reserve);
   std::vector<int> res_end_r1;
   res_end_r1.reserve(memory_reserve);
   std::vector<int> res_start_r2;
   res_start_r2.reserve(memory_reserve);
   std::vector<int> res_end_r2;
   res_end_r2.reserve(memory_reserve);

   std::vector<int> res_r1_avg;
   res_r1_avg.reserve(memory_reserve);
   std::vector<int> res_r2_avg;
   res_r2_avg.reserve(memory_reserve);


   //pileups is a vector of all positions and counts
   int r1_length = end_r1.size();
   //iterate through the vector of end positions
   for (int i = 0; i < r1_length; i++) {
      int r1_start = start_r1[i];
      int r1_end = end_r1[i];
      int r2_start = start_r2[i];
      int r2_end = end_r2[i];
      //set read length to be read end - read start
      //set start to be read_start_vec at pos i
      int read_r1_length = end_r1[i] - start_r1[i] + 1;
      int read_r2_length = end_r2[i] - start_r2[i] + 1;
      //create results vector for indexes to be returned
      std::vector<int> result_vector(read_r1_length);
      int r1_running_total = 0;
      int r2_running_total = 0;


      // Assumes start_r1 and start_r2 are of same length
      // Need to confirm
      for(int j = start_r1[i], k = start_r2[i]; j <= end_r1[i]; j++, k++){
         //create iterator for r1
         std::vector<int>::iterator it_r1;
         std::vector<int>::iterator it_r2;
         it_r1 = std::find(dtpos.begin(), dtpos.end(), j);
         it_r2 = std::find(dtpos.begin(), dtpos.end(), k);
         // Might need to change this to check before incrementing total instead of continuing here
         int index_r1, index_r2;
         if (it_r1 != dtpos.end()) {
            index_r1 = it_r1 - dtpos.begin();
            r1_running_total += dtcount[index_r1];
         }
         if (it_r2 != dtpos.end()) {
            index_r2 = it_r2 - dtpos.begin();
            r2_running_total += dtcount[index_r2];
         }
      }

      int r1_count_average, r2_count_average;
      if (r1_running_total != 0) {
         r1_count_average = r1_running_total / read_r1_length;
         res_start_r1.emplace_back(r1_start);
         res_end_r1.emplace_back(r1_end);
         res_r1_avg.emplace_back(r1_count_average);
      }
      if (r2_running_total != 0) {
         r2_count_average = r2_running_total / read_r2_length;
         res_start_r2.emplace_back(r2_start);
         res_end_r2.emplace_back(r2_end);
         res_r2_avg.emplace_back(r2_count_average);
      }
   }

   DataFrame df = DataFrame::create(Named("r1_start") = res_start_r1, Named("r1_end") = res_end_r1, Named("r1_count_avg") = res_r1_avg,
                                    Named("r2_start") = res_start_r2, Named("r2_end") = res_end_r2, Named("r2_count_avg") = res_r2_avg);
   return df;
}

//[[Rcpp::export]]
std::map<int, int> vectorsToMap(std::vector<int> &k, std::vector<int> &v) {
   // Assumes vectors of same length
   // Only call from getPileupsMap which should only be called with data.frame input
   // std::cout << "!!Assertion upcoming!!" << std::endl;
   // std::cout << "Checking if dtpos.size() and dtcount.size() are equal..." << std::endl;
   // assert(k.size() == v.size());

   std::map<int, int> m;
   std::transform(k.begin(), k.end(), v.begin(), std::inserter(m, m.end()),
                  [](int a, int b) {
                     return std::make_pair(a, b);
                  });
   return m;
}

//[[Rcpp::export]]
DataFrame getPileupsMap(std::vector<int> dtpos, std::vector<int> dtcount, std::vector<int> start_r1, std::vector<int> end_r1,
                     std::vector<int> start_r2, std::vector<int> end_r2, std::vector<int> count){
   //calculates the average pileup for a read using counts at each nucleotide position
   //takes position from read pileups df, count at that position, start pos of candidate read, end pos of candidate read

   int memory_reserve = start_r1.size();

   std::map<int, int> dt = vectorsToMap(dtpos, dtcount);

   std::vector<int> res_start_r1;
   res_start_r1.reserve(memory_reserve);
   std::vector<int> res_end_r1;
   res_end_r1.reserve(memory_reserve);
   std::vector<int> res_start_r2;
   res_start_r2.reserve(memory_reserve);
   std::vector<int> res_end_r2;
   res_end_r2.reserve(memory_reserve);
   std::vector<int> res_r1_avg;
   res_r1_avg.reserve(memory_reserve);
   std::vector<int> res_r2_avg;
   res_r2_avg.reserve(memory_reserve);

   //pileups is a vector of all positions and counts
   int r1_length = end_r1.size();
   //iterate through the vector of end positions
   for (int i = 0; i < r1_length; i++) {
      int r1_start = start_r1[i];
      int r1_end = end_r1[i];
      int r2_start = start_r2[i];
      int r2_end = end_r2[i];
      int duplicate_count = count[i];

      //set read length to be read end - read start
      //set start to be read_start_vec at pos i
      int read_r1_length = r1_end - r1_start + 1;
      int read_r2_length = r2_end - r2_start + 1;

      int r1_running_total = 0;
      int r2_running_total = 0;

      std::map<int, int>::iterator it;

      for (int j = r1_start; j < r1_end; j++) {
         it = dt.find(j);
         if (it != dt.end()) {
            r1_running_total += it->second;
         }
      }
      for (int k = r2_start; k < r2_end; k++) {
         it = dt.find(k);
         if (it != dt.end()) {
            r2_running_total += it->second;
         }
      }

      int r1_count_average, r2_count_average;
      if (r1_running_total != 0) {
         // This was added in to reduce the number of iterations which was approaching 1 billion
         // Used dplyr to group all the duplicate overlap reads and mutate in a count column to keep track of the number of duplicates
         // That way we only have to iterate through each read range once and multiply the results by the number of duplicates
         r1_running_total *= duplicate_count;
         r1_count_average = r1_running_total / read_r1_length;
         res_start_r1.emplace_back(r1_start);
         res_end_r1.emplace_back(r1_end);
         res_r1_avg.emplace_back(r1_count_average);
      }
      if (r2_running_total != 0) {
         r2_running_total *= duplicate_count;
         r2_count_average = r2_running_total / read_r2_length;
         res_start_r2.emplace_back(r2_start);
         res_end_r2.emplace_back(r2_end);
         res_r2_avg.emplace_back(r2_count_average);
      }
   }

   DataFrame df = DataFrame::create(Named("r1_start") = res_start_r1, Named("r1_end") = res_end_r1, Named("r1_count_avg") = res_r1_avg,
                                    Named("r2_start") = res_start_r2, Named("r2_end") = res_end_r2, Named("r2_count_avg") = res_r2_avg);
   return df;
}

//[[Rcpp::export]]
DataFrame group_helix_res(std::vector<int> x, std::vector<int> y) {
   int X_SIZE = x.size();
   int Y_SIZE = y.size();
   if (X_SIZE == 0) {
      std::string msg = "Input cannot be empty";
      DataFrame df = DataFrame::create(Named("Error") = msg);
      return(df);
   }
   if (X_SIZE != Y_SIZE) {
      std::string msg = "Input arguments must be of the same length";
      DataFrame df = DataFrame::create(Named("Error") = msg);
      return(df);
   }

   int x1, x2, y1, y2, previous_x, previous_y;
   std::vector<int> x_start;
   std::vector<int> x_end;
   std::vector<int> y_start;
   std::vector<int> y_end;

   //store first result, keep iterating, then store last result
   for (int i = 0; i < X_SIZE; i++) {
      if (i == 0) { // When starting a new range of numbers
         x1 = x[i];
         y1 = y[i];
      } else {
         if (x[i] - previous_x > 4) {
            x2 = previous_x;
            y2 = previous_y;
            x_start.push_back(x1);
            x_end.push_back(x2);
            y_start.push_back(y1);
            y_end.push_back(y2);
            x1 = x[i];
            y1 = y[i];
         }
      }
      // Store the current x and y for easy access in next iteration
      previous_x = x[i];
      previous_y = y[i];

      if (i == X_SIZE - 1) { // On last iteration of the loop, store the current results
         x2 = x[i];
         y2 = y[i];
         x_start.push_back(x1);
         x_end.push_back(x2);
         y_start.push_back(y1);
         y_end.push_back(y2);
      }
   }
   std::cout << "CPP IS UPDATED" << std::endl;
   //DataFrame result = DataFrame::create(Named("X Start") = x_start, Named("X End") = x_end, Named("Y Start") = y_start, Named("Y End") = y_end);
   DataFrame result = DataFrame::create(Named("X Start") = x_start, Named("X End") = x_end, Named("Y Start") = y_start, Named("Y End") = y_end);
   return(result);
}


//[[Rcpp::export]]
DataFrame calc_overhangs(std::vector<int> r1_start, std::vector<int> r1_end,
                         std::vector<int> r2_start, std::vector<int> r2_width){
   std::vector<int> shift_vec = {-4,-3,-2,-1,0,1,2,3,4};
   int mut_size = r1_start.size();
   int MEMORY_SIZE = 9;

   std::vector<int> proper_count;
   std::vector<int> improper_count;
   proper_count.reserve(MEMORY_SIZE);
   improper_count.reserve(MEMORY_SIZE);
   for(int i = -4; i <= 4; i++){
      int pcount = 0;
      int icount = 0;

      for(int j = 0; j < mut_size; j++){
         int new_r2_start = r1_start[j] + 2 + i;
         int new_r2_end = new_r2_start + r2_width[j] - 1;
         if((new_r2_end - r1_end[j] == 2)){
            pcount += 1;
         } else {
            icount += 1;
         }
      }
      proper_count.emplace_back(pcount);
      improper_count.emplace_back(icount);
   }
   DataFrame overhangs = DataFrame::create(Named("shift") = shift_vec, Named("proper_count") = proper_count,
                                           Named("improper_count") = improper_count);
   return(overhangs);
}

//[[Rcpp::export]]
DataFrame calc_expand_overhangs(std::vector<int> r1_start, std::vector<int> r1_end,
                         std::vector<int> r2_start, std::vector<int> r2_width){
  std::vector<int> shift_vec = {-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8};
  int mut_size = r1_start.size();
  int MEMORY_SIZE = 17;

  std::vector<int> proper_count;
  std::vector<int> improper_count;
  proper_count.reserve(MEMORY_SIZE);
  improper_count.reserve(MEMORY_SIZE);
  for(int i = -8; i <= 8; i++){
    int pcount = 0;
    int icount = 0;

    for(int j = 0; j < mut_size; j++){
      int new_r2_start = r1_start[j] + 2 + i;
      int new_r2_end = new_r2_start + r2_width[j] - 1;
      if((new_r2_end - r1_end[j] == 2)){
        pcount += 1;
      } else {
        icount += 1;
      }
    }
    proper_count.emplace_back(pcount);
    improper_count.emplace_back(icount);
  }
  DataFrame overhangs = DataFrame::create(Named("shift") = shift_vec, Named("proper_count") = proper_count,
                                          Named("improper_count") = improper_count);
  return(overhangs);
}



//[[Rcpp::export]]
int proper_overlap(int widthx, int widthy) {
   // widthx = 15 widthy = 15
   // widthx - 2 = 13 so p_overlap = 15

   // widthx = 18 widthy = 15
   // widthx - 2 = 16 which is greater than 15
   // p_overlap = 15 (if)

   int proper_overlap;

   if (widthx - 2 > widthy) {
      proper_overlap = widthy;
   } else {
      proper_overlap = widthx - 2;
   }
   return proper_overlap;
}


//[[Rcpp::export]]
int overlap_counts(std::vector<int> f_start, int X_SIZE,
                   std::vector<int> r_end, int Y_SIZE,
                   int overlap) {
   // filter read sizes prior to this
   // For each forward read start position, is there a reverse read rp5 position that is fp5+overlap-1?
   // If so then the forward read gets counted and each reverse read gets counted
   int count = 0;

   for (int i = 0; i < X_SIZE; i++) {
      bool i_status = false;
      int calculated_r_pos = f_start[i] + overlap - 1;
      for (int j = 0; j < Y_SIZE; j++) {
         if (calculated_r_pos == r_end[j]) {
            if(!i_status) i_status = true;
            count++;
         }
      }
      if (i_status) count++;
   }

   return count;

}



//[[Rcpp::export]]
NumericMatrix get_si_overlaps(std::vector<int> fdt_start, std::vector<int> fdt_end, std::vector<int> fwidth,
                              std::vector<int> rdt_start, std::vector<int> rdt_end, std::vector<int> rwidth){

   int m_size = 18;
   NumericMatrix result(m_size);
   //to avoid warnings about signed int to unsigned int conversion
   int f_size = int(fdt_start.size());
   int r_size = int(rdt_start.size());

   for(int i = 15; i <= 32; i++) {
      for(int j = 15; j <= 32; j++) {
         int p_overlap = proper_overlap(i,j);
         int current_count = 0;

         if(i == j){
            std::vector<int> fstart_res;
            std::vector<int> fend_res;
            std::vector<int> rstart_res;
            std::vector<int> rend_res;
            //iterate through fdt
            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == i){
                  fstart_res.push_back(fdt_start[k]);
                  fend_res.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == j){
                  rstart_res.push_back(rdt_start[l]);
                  rend_res.push_back(rdt_end[l]);
               }
            }

            current_count = overlap_counts(fstart_res, fstart_res.size(), rstart_res, rstart_res.size(), p_overlap);
         } else {
            std::vector<int> fstart_res1;
            std::vector<int> fend_res1;
            std::vector<int> rstart_res1;
            std::vector<int> rend_res1;

            std::vector<int> fstart_res2;
            std::vector<int> fend_res2;
            std::vector<int> rstart_res2;
            std::vector<int> rend_res2;

            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == i){
                  fstart_res1.push_back(fdt_start[k]);
                  fend_res1.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == j){
                  rstart_res1.push_back(rdt_start[l]);
                  rend_res1.push_back(rdt_end[l]);
               }
            }
            int c1 = overlap_counts(fstart_res1, fstart_res1.size(), rstart_res1, rstart_res1.size(), p_overlap);

            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == j){
                  fstart_res2.push_back(fdt_start[k]);
                  fend_res2.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == i){
                  rstart_res2.push_back(rdt_start[l]);
                  rend_res2.push_back(rdt_end[l]);
               }
            }
            int c2 = overlap_counts(fstart_res2, fstart_res2.size(), rstart_res2, rstart_res2.size(), p_overlap);

            current_count = c1 + c2;


         }
         result(i - 15,j - 15) = current_count;
      }

   }
   return(result);
}


//[[Rcpp::export]]
DataFrame get_phased_dist(IntegerVector startX, IntegerVector endX,
                          IntegerVector startY, IntegerVector endY, int dfsize){
   // Create an initial memory block size for the vectors
   int memory_reserve = dfsize;

   // Create vectors to hold each new data frame column and reserve initial memory for them
   std::vector<int> res_startx;
   res_startx.reserve(memory_reserve);
   std::vector<int> res_starty;
   res_starty.reserve(memory_reserve);
   std::vector<int> res_dist;
   res_dist.reserve(memory_reserve);
   std::vector<int> res_width_x;
   res_width_x.reserve(memory_reserve);
   std::vector<int> res_width_y;
   res_width_y.reserve(memory_reserve);

   // Loop over the input in a for loop using dfsize as the limit
   for (int i = 0; i < dfsize; i++){
      int start_x = startX[i];
      int end_x = endX[i];
      int w1 = endX[i] - startX[i] + 1;
      int start_y = startY[i];

      int w2 = endY[i] - startY[i] + 1;
      // compute the distance between the values
      int dist = end_x - start_y;
      // matched sequences cannot be overlapping
      if ((start_x < start_y) && (end_x < start_y)){
         res_startx.emplace_back(start_x);
         res_starty.emplace_back(start_y);
         res_dist.emplace_back(dist);
         res_width_x.emplace_back(w1);
         res_width_y.emplace_back(w2);
      }

   }


   // After loops create a new dataframe and assign each column to one of the vectors we created
   DataFrame df = DataFrame::create(Named("start_x") = res_startx,
                                    Named("start_y") = res_starty, Named("width_x") = res_width_x,
                                    Named("width_y") = res_width_y, Named("dist") = res_dist);
   return df;

}




//[[Rcpp::export]]

DataFrame make_count_table(std::vector<int> fdt_start, std::vector<int> fdt_end, std::vector<int> fwidth,
                           std::vector<int> rdt_start, std::vector<int> rdt_end, std::vector<int> rwidth){
//function to calculate counts of each overlap to be used for Z score
//calls the Rcpp function overlap_counts
   int res_size = 27;

   std::vector<int> overlap_res(res_size);
   std::vector<int> counts_res(res_size);

   int f_size = int(fdt_start.size());
   int r_size = int(rdt_start.size());
   for(int maxquery = 18; maxquery <= 30; maxquery ++){
      for(int overlap = 4; overlap <= 30; overlap ++){

         std::vector<int> fstart_res;
         std::vector<int> fend_res;
         std::vector<int> rstart_res;
         std::vector<int> rend_res;

         for(int k = 0; k < f_size; k++){
            //get reads of size i
            if(fwidth[k] == maxquery){
               fstart_res.push_back(fdt_start[k]);
               fend_res.push_back(fdt_end[k]);
            }
         }
         for(int l = 0; l < r_size; l++){
            if(rwidth[l] == maxquery){
               rstart_res.push_back(rdt_start[l]);
               rend_res.push_back(rdt_end[l]);
            }
         }
         overlap_res[overlap - 4] = overlap;
         int c1 = overlap_counts(fstart_res, fstart_res.size(), rend_res, rend_res.size(), overlap);
         counts_res[overlap - 4] = counts_res[overlap - 4] + c1;
      }

   }
   DataFrame df = DataFrame::create(Named("overlap") = overlap_res, Named("count")= counts_res);
   return(df);
}




//[[Rcpp::export]]
NumericMatrix get_pi_overlaps(std::vector<int> fdt_start, std::vector<int> fdt_end, std::vector<int> fwidth,
                              std::vector<int> rdt_end, std::vector<int> rdt_start, std::vector<int> rwidth){

   int m_size = 18;
   NumericMatrix result(m_size);
   //to avoid warnings about signed int to unsigned int conversion
   int f_size = int(fdt_start.size());
   int r_size = int(rdt_start.size());

   for(int i = 15; i <= 32; i++) {
      for(int j = 15; j <= 32; j++) {
         int p_overlap = 10;
         int current_count = 0;

         if(i == j){
            std::vector<int> fstart_res;
            std::vector<int> fend_res;
            std::vector<int> rstart_res;
            std::vector<int> rend_res;
            //iterate through fdt
            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == i){
                  fstart_res.push_back(fdt_start[k]);
                  fend_res.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == j){
                  rstart_res.push_back(rdt_start[l]);
                  rend_res.push_back(rdt_end[l]);
               }
            }
            current_count = overlap_counts(fstart_res, fstart_res.size(), rend_res, rend_res.size(), p_overlap);

         } else {
            std::vector<int> fstart_res1;
            std::vector<int> fend_res1;
            std::vector<int> rstart_res1;
            std::vector<int> rend_res1;

            std::vector<int> fstart_res2;
            std::vector<int> fend_res2;
            std::vector<int> rstart_res2;
            std::vector<int> rend_res2;

            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == i){
                  fstart_res1.push_back(fdt_start[k]);
                  fend_res1.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == j){
                  rstart_res1.push_back(rdt_start[l]);
                  rend_res1.push_back(rdt_end[l]);
               }
            }
            int c1 = overlap_counts(fstart_res1, fstart_res1.size(), rend_res1, rend_res1.size(), p_overlap);


            for(int k = 0; k < f_size; k++){
               //get reads of size i
               if(fwidth[k] == j){
                  fstart_res2.push_back(fdt_start[k]);
                  fend_res2.push_back(fdt_end[k]);
               }
            }
            for(int l = 0; l < r_size; l++){
               if(rwidth[l] == i){
                  rstart_res2.push_back(rdt_start[l]);
                  rend_res2.push_back(rdt_end[l]);
               }
            }
            int c2 = overlap_counts(fstart_res2, fstart_res2.size(), rend_res2, rend_res2.size(), p_overlap);

            current_count = c1 + c2;

         }
         result(i - 15,j - 15) = current_count;
      }

   }
   return(result);
}


//[[Rcpp::export]]
DataFrame getLoopPileupsCPP(std::vector<int> r1Start, std::vector<int> r1Stop,
                            std::vector<int> lStart, std::vector<int> lStop,
                            std::vector<int> r2Start, std::vector<int> r2Stop,
                            std::vector<int> dt_col_1, std::vector<int> dt_col_2,
                            int total_count) {
  // Expects 9 vectors of integers. The first 6 vectors and the 7th and 8th vectors should have the same number of elements.
  unsigned int v1_length = r1Start.size();
  unsigned int dt_length = dt_col_1.size();

  // Confirm first 6 vectors have same number of elements in case this ever gets manually called
  if (v1_length != r1Stop.size()  ||
      v1_length != lStart.size()  ||
      v1_length != lStop.size()   ||
      v1_length != r2Start.size() ||
      v1_length != r2Stop.size())   {
    std::string msg = "Error, all inputs must have the same number of elements. Consider calling from a data.frame.";
    DataFrame df = DataFrame::create(Named("Error") = msg);
    return(df);
  }

  // Confirm the 7th and 8th vectors have the same number of elements
  if (dt_length != dt_col_2.size()) {
    std::string msg = "Error, dt_col_1 and dt_col_2 should have the same number of elements. Consider calling from a data.frame.";
    DataFrame df = DataFrame::create(Named("Error") = msg);
  }

  std::vector<int> r1_start;
  std::vector<int> r1_stop;
  std::vector<int> l_start;
  std::vector<int> l_stop;
  std::vector<int> r2_start;
  std::vector<int> r2_stop;
  r1_start.reserve(v1_length);
  r1_stop.reserve(v1_length);
  l_start.reserve(v1_length);
  l_stop.reserve(v1_length);
  r2_start.reserve(v1_length);
  r2_stop.reserve(v1_length);

  for (unsigned int i = 0; i < v1_length; i++) {
    int cur_l_start = lStart[i];
    int cur_l_stop = lStop[i];
    int count = 0;

    // Check to see if dt_col_1 is in range of the loop start and stop
    // If so, add the integer in dt_col_2 to count
    for (unsigned int j = 0; j < dt_length; j++) {
      if (dt_col_1[j] >= cur_l_start &&
          dt_col_1[j] <= cur_l_stop) {
        count += dt_col_2[j];
      }
    }

    //total_count += count;
    if (count/(total_count + count) < 0.05) {
      r1_start.emplace_back(r1Start[i]);
      r1_stop.emplace_back(r1Stop[i]);
      l_start.emplace_back(cur_l_start);
      l_stop.emplace_back(cur_l_stop);
      r2_start.emplace_back(r2Start[i]);
      r2_stop.emplace_back(r2Stop[i]);
    }
  }

  DataFrame df = DataFrame::create(Named("r1_start") = r1_start, Named("r1_stop") = r1_stop,
                                   Named("lstart") = l_start, Named("lstop") = l_stop,
                                   Named("r2_start") = r2_start, Named("r2_stop") = r2_stop);
  return df;
}

