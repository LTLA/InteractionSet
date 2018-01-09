#include "overlap_utils.h"

/*****************************************************************************************************
 * Another function that links two regions together. This does something quite similar
 * to detect_paired_olaps, but whereas that function requires both anchor regions to
 * overlap the regions from the same subject pair, here we only require that the
 * anchor regions overlap one region from either set. We then store all possible 
 * combinations of regions that are mediated by our interactions.
 *****************************************************************************************************/

SEXP expand_pair_links(SEXP anchor1, SEXP anchor2, SEXP querystarts1, SEXP queryends1, SEXP subject1, SEXP nsubjects1, 
        SEXP querystarts2, SEXP queryends2, SEXP subject2, SEXP nsubjects2, SEXP sameness, SEXP use_both) {
    BEGIN_RCPP

    Rcpp::IntegerVector _a1(anchor1), _a2(anchor2);
    const int Npairs = _a1.size();
    if (Npairs != _a2.size()) { 
        throw std::runtime_error("anchor vectors must be of equal length");
    }
   
    // Checking start/end/subject indices for region 1.
    Rcpp::IntegerVector _qs1(querystarts1), _qe1(queryends1), _subj1(subject1), _nsubjects1(nsubjects1);
    if (_nsubjects1.size()!=1) { 
        throw std::runtime_error("total number of subjects (1) must be an integer scalar");
    }
    const int Ns_all1 = _nsubjects1[0];
    check_indices(_qs1, _qe1, _subj1, Ns_all1);

    // Checking start/end/subject indices for region 2.
    Rcpp::IntegerVector _qs2(querystarts2), _qe2(queryends2), _subj2(subject2), _nsubjects2(nsubjects2);
    if (_nsubjects2.size()!=1) { 
        throw std::runtime_error("total number of subjects (2) must be an integer scalar");
    }
    const int Ns_all2 = _nsubjects2[0];
    check_indices(_qs2, _qe2, _subj2, Ns_all2);

    const int Nq=_qs1.size();
    if (Nq!=_qs2.size()) { 
        throw std::runtime_error("start/end vectors should be the same length for regions 1 and 2");
    }

    // Checking other parameters.
    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);           

    Rcpp::LogicalVector _sameness(sameness);
    if (_sameness.size()!=1) { 
        throw std::runtime_error("same region specification should be a logical scalar"); 
    }
    const bool is_same=_sameness[0];
    if (is_same) { // No point examining the flipped ones.
        true_mode_end=true_mode_start+1; 
    } 
  
    // Setting up the intermediate structures.
    typedef std::pair<int, int> link;
    std::set<link> currently_active;
    std::set<link>::const_iterator itca;
    std::deque<link> stored_inactive;
    std::deque<int> interactions;

    // Running through each pair and seeing what it links together.
    for (int curpair=0; curpair<Npairs; ++curpair) {
        const int& a1=_a1[curpair];
        const int& a2=_a2[curpair];
        const int& maxmode = (a1==a2 ? true_mode_start+1 : true_mode_end);

        // Repeating with switched anchors, if the query sets are not the same.
        for (int mode=true_mode_start; mode<maxmode; ++mode) { 
            int curq1=0, curq2=0;
            if (mode==0) { 
                curq1 = a1;
                curq2 = a2;
                if (curq1 >= Nq || curq1 < 0 || curq1==NA_INTEGER) { 
                    throw std::runtime_error("region index (1) out of bounds"); 
                }
                if (curq2 >= Nq || curq2 < 0 || curq2==NA_INTEGER) { 
                    throw std::runtime_error("region index (2) out of bounds"); 
                }
            } else {
                curq2 = a1;
                curq1 = a2;
            }

            /* Storing all combinations associated with this pair. Avoiding redundant combinations
             * for self-linking (we can't use a simple rule like curindex2 < curindex1 to restrict 
             * the loop, as the second anchor can overlap regions above the first anchor if the 
             * latter is nested within the former).
             */
            const int& qs1=_qs1[curq1];
            const int& qe1=_qe1[curq1];
            const int& qs2=_qs2[curq2];
            const int& qe2=_qe2[curq2];
            for (int q1=qs1; q1<qe1; ++q1) {
                const int& s1=_subj1[q1];
                for (int q2=qs2; q2<qe2; ++q2) { 
                    const int& s2=_subj2[q2];

                    if (is_same && s1 < s2) {
                        currently_active.insert(link(s2, s1));
                    } else {
                        currently_active.insert(link(s1, s2));
                    }
                }
            }
        }

        // Relieving the set by storing all inactive entries (automatically sorted as well).
        stored_inactive.insert(stored_inactive.end(), currently_active.begin(), currently_active.end());
        interactions.resize(interactions.size() + currently_active.size(), curpair);
        currently_active.clear();
    }

    // Popping back a list of information.
    Rcpp::IntegerVector out_inter(interactions.begin(), interactions.end());
    const int total_entries=interactions.size();
    Rcpp::IntegerVector out_first(total_entries), out_second(total_entries);
    Rcpp::IntegerVector::iterator ofIt=out_first.begin(), osIt=out_second.begin();

    for (int curdex=0; curdex < total_entries; ++curdex, ++ofIt, ++osIt) {
        *ofIt=stored_inactive[curdex].first;
        *osIt=stored_inactive[curdex].second;
    }
    
    return Rcpp::List::create(out_inter, out_first, out_second);
    END_RCPP
}

