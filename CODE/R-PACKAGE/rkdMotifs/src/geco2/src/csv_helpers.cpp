#include <csv_iterator.hpp>

void csv::details::checkIteratorRange(const strIt& curr, const strIt& end){
    if(curr == end){
        throw(std::out_of_range("Not enough parameter building the tuple"));
    }
}
