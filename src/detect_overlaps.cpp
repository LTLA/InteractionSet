#include "overlap_utils.h"

/*****************************************************************************************************
 * Useful classes to get different outputs from the same code. *
 *****************************************************************************************************/

class output_store {
public:
    output_store() {}
    virtual ~output_store() {}
    virtual void prime(int, int)=0;
    virtual void acknowledge(int, int)=0;
    virtual void postprocess()=0;
    virtual Rcpp::RObject generate() const=0;
    virtual bool quit () const { 
        return false; 
    }
};

class expanded_overlap : public output_store { // Stores all new query:subject relations.
public:
    expanded_overlap() : just_added(0) {};
    ~expanded_overlap() {};
    void prime(int nq, int ns) {
        (void)nq; (void)ns;
    }
    void acknowledge(int q, int s) {
        new_query.push_back(q);
        new_subject.push_back(s);
        ++just_added;
    }
    void postprocess() { // Sorting for output.
        std::sort(new_subject.end()-just_added, new_subject.end());
        just_added=0;
    }
    Rcpp::RObject generate () const {
        return Rcpp::List::create(
            Rcpp::IntegerVector(new_query.begin(), new_query.end()),
            Rcpp::IntegerVector(new_subject.begin(), new_subject.end())
        );
    }
private:
    std::deque<int> new_query, new_subject;
    int just_added;
};

class single_query_overlap : public output_store { // Stores a 1:1 query->subject relationship.
public:
    single_query_overlap() {}
    ~single_query_overlap() {}
    void prime(int nq, int ns) {
        (void)ns;
        nquery=nq;
        tosubject.resize(nq, NA_INTEGER);    
    }
    void acknowledge(int q, int s) {
        if (q >= nquery) { 
            throw std::runtime_error("requested query index out of range"); 
        }
        tosubject[q]=s;
    }
    void postprocess() {}
    Rcpp::RObject generate () const {
        return Rcpp::IntegerVector(tosubject.begin(), tosubject.end());
    } 
    bool quit () const { 
        return false; 
    }  
protected:
    int nquery;
    std::deque<int> tosubject;
};

class first_query_overlap : public single_query_overlap { // Stores the first.
public:
    first_query_overlap() {}
    ~first_query_overlap() {}
    void acknowledge(int q, int s) {
        if (q >= nquery) { 
            throw std::runtime_error("requested query index out of range"); 
        }
        if (tosubject[q]==NA_INTEGER || tosubject[q] > s) {
            tosubject[q]=s;
        }
    }
};

class last_query_overlap : public single_query_overlap { // Stores the last.
public:
    last_query_overlap() {}
    ~last_query_overlap() {}
    void acknowledge(int q, int s) {
        if (q >= nquery) { 
            throw std::runtime_error("requested query index out of range"); 
        }
        if (tosubject[q]==NA_INTEGER || tosubject[q] < s) {
            tosubject[q]=s;
        }
    }
};

class arbitrary_query_overlap : public single_query_overlap { // Stores any (and quits once it finds one).
    public:    
    arbitrary_query_overlap() {}
    ~arbitrary_query_overlap() {}
    bool quit () const { 
        return true; 
    } 
};

class single_subject_overlap : public output_store { // Stores a 1:1 subject->query relationship.
public:
    single_subject_overlap() {}
    ~single_subject_overlap() {}
    void prime(int nq, int ns) {
        (void)nq;
        nsubject=ns;
        toquery.resize(ns, NA_INTEGER);    
    };
    void acknowledge(int q, int s) {
        if (s >= nsubject) { 
            throw std::runtime_error("requested subject index out of range"); 
        }
        toquery[s]=q;
    }
    void postprocess() {}
    Rcpp::RObject generate () const {
        return Rcpp::IntegerVector(toquery.begin(), toquery.end());
    } 
    bool quit () const { 
        return false; 
    } 
protected:
    int nsubject;
    std::deque<int> toquery;
};

class first_subject_overlap : public single_subject_overlap { // Stores the first.
public:
    first_subject_overlap() {}
    ~first_subject_overlap() {}
    void acknowledge(int q, int s) {
        if (s >= nsubject) { 
            throw std::runtime_error("requested subject index out of range"); 
        }
        if (toquery[s]==NA_INTEGER || toquery[s] > q) {
            toquery[s]=q;
        }
    }
};

class last_subject_overlap : public single_subject_overlap { // Stores the last.
public:
    last_subject_overlap() {}
    ~last_subject_overlap() {}
    void acknowledge(int q, int s) {
        if (s >= nsubject) { 
            throw std::runtime_error("requested subject index out of range"); 
        }
        if (toquery[s]==NA_INTEGER || toquery[s] < q) {
            toquery[s]=q;
        }
    }
};

class query_count_overlap : public output_store { // Stores number of times 'query' was hit.
public:
    query_count_overlap() : nquery(0) {}
    ~query_count_overlap() {}
    void prime(int nq, int ns) {
        (void)ns;
        nquery=nq;
        query_hit.resize(nq);    
    }
    void acknowledge(int q, int s) {
        (void)s;
        if (q >= nquery) { 
            throw std::runtime_error("requested query index out of range"); 
        }
        ++query_hit[q];
    }
    void postprocess() {}
    Rcpp::RObject generate () const {
        return Rcpp::IntegerVector(query_hit.begin(), query_hit.end());
    } 
private:
    int nquery;
    std::deque<int> query_hit;
};

class subject_count_overlap : public output_store { // Stores number of times 'subject' was hit.
public:
    subject_count_overlap() : nsubject(0) {}
    ~subject_count_overlap() {}
    void prime(int nq, int ns) {
        (void)nq;
        nsubject=ns;
        subject_hit.resize(ns);    
    }
    void acknowledge(int q, int s) {
        (void)q;
        if (s >= nsubject) { 
            throw std::runtime_error("requested subject index out of range"); 
        }
        ++subject_hit[s];
    }
    void postprocess() {}
    Rcpp::RObject generate () const {
        return Rcpp::IntegerVector(subject_hit.begin(), subject_hit.end());
    } 
private:
    int nsubject;
    std::deque<int> subject_hit;
};

/*****************************************************************************************************
 * Base function, to detect all overlaps between linear ranges and either interacting loci in a pair *
 *****************************************************************************************************/

void help_add_current_overlaps(const int& true_mode_start, const int& true_mode_end,
        const int& curpair, 
        const Rcpp::IntegerVector& anchor1, const Rcpp::IntegerVector& anchor2,
        const Rcpp::IntegerVector& query_start, const Rcpp::IntegerVector& query_end,
        const Rcpp::IntegerVector& subject_index, 
        output_store* output, std::vector<int>& latest_pair) {
    /* This helper function checks all overlaps involving the interaction indexed by 'curpair'.
     * I factorized it out of the main function below to avoid using 'goto's upon the quit() call.
     */

    const int& a1=anchor1[curpair];
    const int& a2=anchor2[curpair];
    const int maxmode=(a1==a2 ? true_mode_start+1 : true_mode_end); // just check one or the other, if they're the same.
    const int Nq=query_start.size();

    for (int mode=true_mode_start; mode<maxmode; ++mode) { 
        const int& curq = (mode == 0 ? a1 : a2);
        if (curq >= Nq || curq < 0 || curq==NA_INTEGER) { 
            throw std::runtime_error("region index out of bounds"); 
        }

        const int& qs=query_start[curq];
        const int& qe=query_end[curq];
        for (int q=qs; q<qe; ++q) {
            const int& curs=subject_index[q];
            if (latest_pair[curs] < curpair) { 
                output->acknowledge(curpair, curs);
                latest_pair[curs] = curpair;

                if (output->quit()) { // If we just want any hit, we go to the next 'curpair'.
                    return;
                }
            }
        }
    }
    return;
}

void detect_olaps(output_store* output, SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both) {
    const Rcpp::IntegerVector a1(anchor1), a2(anchor2);
    const int Npairs = a1.size();
    if (Npairs != a2.size()) { 
        throw std::runtime_error("anchor vectors must be of equal length"); 
    } 
   
    const Rcpp::IntegerVector qs(querystarts), qe(queryends), subj(subject), _nsubjects(nsubjects);
    if (_nsubjects.size()!=1) { 
        throw std::runtime_error("total number of subjects must be an integer scalar");
    }
    const int Ns_all = _nsubjects[0];
    check_indices(qs, qe, subj, Ns_all);
    
    // Checking which regions to use for the overlap.
    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);

    // Detecting overlaps between each pair of anchor regions.
    output->prime(Npairs, Ns_all);
    std::vector<int> latest_pair(Ns_all, -1);

    for (int curpair=0; curpair<Npairs; ++curpair) {
        help_add_current_overlaps(true_mode_start, true_mode_end, 
                curpair, a1, a2, qs, qe, subj, output, latest_pair);
        output->postprocess();
    }
    return;
}

/*****************************************************************************************************
 * Base function, to detect all 2D overlaps between two GInteractions objects.
 * This is pretty complicated:
 *
 *   anchor1, anchor2 hold the anchor indices of the query GInteractions.
 *   querystarts, queryend hold the first and one-after-last row of same-query runs in the Hits object obtained by findOverlaps(regions(query), regions(subject))
 *   subject holds the subjectHits of the aforementioned Hits object.
 *   next_anchor_start1, next_anchor_end1 holds the first and one-after-last position of same-value runs in the sorted subject@anchor1.
 *   next_id1 holds the order of subject@anchor1.
 *   next_anchor_start2, next_anchor_end2 holds the first and one-after-last position of same value runs in the sorted subject@anchor2.
 *   next_id2 holds the order of subject@anchor2.
 *   next_num_pairs holds the total number of pairs in the subject.
 * 
 * The idea is to:
 *   1) go through each query pair;
 *   2) retrieve the regions(subject) overlapping each anchor query;
 *   3) identify all subject anchors matching each retrieved regions(subject)
 *   4) cross-reference all identified subject anchors 1 and 2 to identify 2D overlaps.
 *****************************************************************************************************/

void help_add_current_paired_overlaps(const int& true_mode_start, const int& true_mode_end,
        const int& curpair, 
        const Rcpp::IntegerVector& anchor1, const Rcpp::IntegerVector& anchor2,
        const Rcpp::IntegerVector& query_start, const Rcpp::IntegerVector& query_end,
        const Rcpp::IntegerVector& next_anchor_start1, const Rcpp::IntegerVector& next_anchor_end1, const Rcpp::IntegerVector& next_id1,
        const Rcpp::IntegerVector& next_anchor_start2, const Rcpp::IntegerVector& next_anchor_end2, const Rcpp::IntegerVector& next_id2,
        const Rcpp::IntegerVector& subject_index, 
        output_store* output, 
        int* latest_pair_A, int* is_complete_A,
        int* latest_pair_B, int* is_complete_B) {
    /* This helper function checks all overlaps involving the interaction indexed by 'curpair'.
     * Factorized out of the main function below to avoid using 'goto's upon the quit() call.
     */

    const int& a1=anchor1[curpair];
    const int& a2=anchor2[curpair];
    const int maxmode = (a1==a2 ? true_mode_start+1 : true_mode_end);
    const int Nq=query_start.size();

    /* Checking whether the first and second anchor overlaps anything in the opposing query sets.
     * Doing this twice; first and second anchors to the first and second query sets (A, mode=0), then
     * the first and second anchors to the second and first query sets (B, mode=1).
     */
    for (int mode=true_mode_start; mode<maxmode; ++mode) { 
        int curq1=0, curq2=0;
        int* latest_pair=NULL;
        int* is_complete=NULL;

        if (mode==0) { 
            curq1 = a1;
            curq2 = a2;
            if (curq1 >= Nq || curq1 < 0 || curq1==NA_INTEGER) { 
                throw std::runtime_error("region index (1) out of bounds"); 
            }
            if (curq2 >= Nq || curq2 < 0 || curq2==NA_INTEGER) { 
                throw std::runtime_error("region index (2) out of bounds"); 
            }
            latest_pair = latest_pair_A;
            is_complete = is_complete_A;
        } else {
            curq2 = a1;
            curq1 = a2;
            latest_pair = latest_pair_B;
            is_complete = is_complete_B;
        }

        const int& qs1=query_start[curq1];
        const int& qe1=query_end[curq1];
        for (int q=qs1; q<qe1; ++q) {
            const int& curs=subject_index[q];
            const int& nas=next_anchor_start1[curs];
            const int& nae=next_anchor_end1[curs];

            for (int n=nas; n<nae; ++n) {
                const int& curid=next_id1[n];

                /* An overlap with element 'cur_nextid' has already been added from the "A" cycle, if 
                 * is_complete[cur_nextid]=true and latest_pair[cur_nextid] is at the current query index.
                 * Otherwise, updating latest_pair if it hasn't already been updated - setting 
                 * latest_pair to the query index to indicate that the first anchor region is overlapped,
                 * but also setting is_complete to false to indicate that the overlap is not complete.
                 */
                if (mode!=0 && latest_pair_A[curid] == curpair && is_complete_A[curid]) { 
                    continue; 
                } 
                if (latest_pair[curid] < curpair) { 
                    latest_pair[curid] = curpair;
                    is_complete[curid] = false;
                }
            }
        }

        const int& qs2=query_start[curq2];
        const int& qe2=query_end[curq2];
        for (int q=qs2; q<qe2; ++q) {
            const int& curs=subject_index[q];
            const int& nas=next_anchor_start2[curs];
            const int& nae=next_anchor_end2[curs];

            for (int n=nas; n<nae; ++n) {
                const int& curid=next_id2[n];

                /* Again, checking if overlap has already been added from the "A" cycle. Otherwise, we only add
                 * an overlap if anchor region 1 was overlapped by 'cur_nextid', as indicated by latest_pair
                 * (and is_complete is false, to avoid re-adding something that was added in this cycle).
                 */
                if (mode!=0 && latest_pair_A[curid] == curpair && is_complete_A[curid]) { 
                    continue; 
                }
                if (latest_pair[curid] == curpair && !is_complete[curid]) {
                    output->acknowledge(curpair, curid);
                    is_complete[curid] = true;
            
                    if (output->quit()) { // If we just want any hit, we go to the next 'curpair'.
                        return;
                    }
                }
            }
        }
    }
    return;
}
 
void detect_paired_olaps(output_store* output, SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject, 
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP use_both) {

    const Rcpp::IntegerVector a1(anchor1), a2(anchor2);
    const int Npairs = a1.size();
    if (Npairs != a2.size()) {
        throw std::runtime_error("anchor vectors must be of equal length"); 
    } 

    const Rcpp::IntegerVector qs(querystarts), qe(queryends), subj(subject);
    const Rcpp::IntegerVector nas1(next_anchor_start1), nae1(next_anchor_end1), nid1(next_id1);
    const Rcpp::IntegerVector nas2(next_anchor_start2), nae2(next_anchor_end2), nid2(next_id2);

    const int Nnp = nid1.size();
    if (Nnp!=nid2.size()) { 
        throw std::runtime_error("number of next IDs is not equal to specified number of pairs"); 
    }

    // Check indices.
    const int Nas=nas1.size();
    if (Nas != nas2.size()) { 
        throw std::runtime_error("run vectors must be of the same length for anchors 1 and 2");
    }
    check_indices(qs, qe, subject, Nas);
    check_indices(nas1, nae1, nid1, Nnp);
    check_indices(nas2, nae2, nid2, Nnp);

    // Check mode.
    int true_mode_start, true_mode_end;
    set_mode_values(use_both, true_mode_start, true_mode_end);

    // Setting up logging arrays. 
    output->prime(Npairs, Nnp);
    std::vector<int> latest_pair_A(Nnp, -1), latest_pair_B(Nnp, -1),
        is_complete_A(Nnp, true), is_complete_B(Nnp, true);

    for (int curpair=0; curpair<Npairs; ++curpair) {
        help_add_current_paired_overlaps(true_mode_start, true_mode_end, 
                curpair, a1, a2, qs, qe, nas1, nae1, nid1, nas2, nae2, nid2,
                subj, output,
                latest_pair_A.data(), is_complete_A.data(),
                latest_pair_B.data(), is_complete_B.data());
        output->postprocess();
    }

    return;
}

/*****************************************************************************************************
 * Functions that actually get called from R.
 *****************************************************************************************************/

std::unique_ptr<output_store> choose_output_type(SEXP select, SEXP GIquery) {
    const Rcpp::StringVector _select(select);
    if (_select.size()!=1) { 
        throw std::runtime_error("'select' specifier should be a single string"); 
    }
    const Rcpp::String selection(_select[0]);
    const char* selstring=selection.get_cstring();
    
    const Rcpp::LogicalVector _GIquery(GIquery);
    if (_GIquery.size()!=1) { 
        throw std::runtime_error("'GIquery' specifier should be a logical scalar"); 
    }
    const bool giq=_GIquery[0];

    if (std::strcmp(selstring, "all")==0) {
        return std::unique_ptr<output_store>(new expanded_overlap);
    } else if (std::strcmp(selstring, "first")==0) {
        if (giq) {
            return std::unique_ptr<output_store>(new first_query_overlap);
        } else {
            return std::unique_ptr<output_store>(new first_subject_overlap);
        }
    } else if (std::strcmp(selstring, "last")==0) {
        if (giq) {
            return std::unique_ptr<output_store>(new last_query_overlap);
        } else {
            return std::unique_ptr<output_store>(new last_subject_overlap);
        }
    } else if (std::strcmp(selstring, "arbitrary")==0) {
        if (giq) {
            return std::unique_ptr<output_store>(new arbitrary_query_overlap);
        } else {
            /* Unfortunately, this CANNOT be sped up via quit(), because the loop is done with 
             * respect to the left GInteractions-as-query. Quitting would only record an arbitrary
             * overlap with respect to the query, not with respect to the subject. Thus, some 
             * subjects that might have been overlapped will be missed when you quit.
             */
            return std::unique_ptr<output_store>(new first_subject_overlap); 
        }
    } else if (std::strcmp(selstring, "count")==0) {
        if (giq) {
            return std::unique_ptr<output_store>(new query_count_overlap);
        } else {
            return std::unique_ptr<output_store>(new subject_count_overlap);
        }
    } 
    throw std::runtime_error("'select' should be 'all', 'first', 'last', 'arbitrary', or 'count'");
}

SEXP linear_olaps(SEXP anchor1, SEXP anchor2, SEXP querystarts, SEXP queryends, SEXP subject, SEXP nsubjects, SEXP use_both, SEXP select, SEXP GIquery) {
    BEGIN_RCPP
    auto x=choose_output_type(select, GIquery);
    detect_olaps(x.get(), anchor1, anchor2, querystarts, queryends, subject, nsubjects, use_both);
    return x->generate();
    END_RCPP
}

SEXP paired_olaps(SEXP anchor1, SEXP anchor2, 
        SEXP querystarts, SEXP queryends, SEXP subject,
        SEXP next_anchor_start1, SEXP next_anchor_end1, SEXP next_id1,
        SEXP next_anchor_start2, SEXP next_anchor_end2, SEXP next_id2,
        SEXP use_both, SEXP select) {
    BEGIN_RCPP
    auto x=choose_output_type(select, Rf_ScalarLogical(1));
    detect_paired_olaps(x.get(), anchor1, anchor2, 
        querystarts, queryends, subject,
        next_anchor_start1, next_anchor_end1, next_id1,
        next_anchor_start2, next_anchor_end2, next_id2,
        use_both);
    return x->generate();
    END_RCPP
}


