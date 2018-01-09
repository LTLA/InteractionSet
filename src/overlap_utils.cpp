#include "overlap_utils.h"

/* Checks the indices specifying runs of the same query in a SortedHits,
 * as well as the subject indices at the intervals spanned by those runs.
 */

void check_indices (const Rcpp::IntegerVector& query_start, const Rcpp::IntegerVector& query_end, 
                    const Rcpp::IntegerVector& subject_index, int Nsubjects) {

    const int Nq=query_start.size();
    if (Nq!=query_end.size()) {
        throw std::runtime_error("vectors of run starts/ends of must have the same length");
    }
    const int Ns=subject_index.size();
    
    // Checking query inputs.
    auto qsIt=query_start.begin(), qeIt=query_end.begin();
    for (int checkdex=0; checkdex < Nq; ++checkdex, ++qsIt, ++qeIt) { 
        const int& curqs =*qsIt;
        const int& curqe =*qeIt;
        if (curqs==NA_INTEGER || curqe==NA_INTEGER) { 
            throw std::runtime_error("indices must be finite integers"); 
        }
        if (curqs >= curqe)  { // Empty range, no need to check actual values.
            continue; 
        } 
        if (curqs >= Ns || curqs < 0) { 
            throw std::runtime_error("start index out of bounds"); 
        }
        if (curqe > Ns || curqe < 0) { 
            throw std::runtime_error("end index out of bounds"); 
        }
    }
    
    // Checking subject inputs.
    if (Nsubjects < 0) { 
        throw std::runtime_error("total number of subject indices must be non-negative"); 
    }
    for (const auto& curs : subject_index) { 
        if (curs >= Nsubjects || curs < 0 || curs==NA_INTEGER) { 
            throw std::runtime_error("subject index out of bounds"); 
        }
    }
    return;
}

/* Sets the mode values - whether or not to check both orientations,
 * just the first->first/second->second (iteration 0) or just the
 * first->second/second->first (iteration 1).
 */

void set_mode_values (Rcpp::IntegerVector use_both, int& start, int& end) {
    if (use_both.size()!=1) { 
        throw std::runtime_error("'use_both' specifier should be an integer scalar"); 
    }
    switch (use_both[0]) { 
        case 1:
            start=0; end=2; break;
        case 2:
            start=0; end=1; break;
        case 3:
            start=1; end=2; break;
        default:
            throw std::runtime_error("invalid specification for 'use_both'");
    }
    return;
}  
