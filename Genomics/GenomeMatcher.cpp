#include "provided.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>

#include "Trie.h"
using namespace std;

bool comparePairByGenome(pair<int, int> p1, pair<int, int> p2){
    if(p1.first == p2.first){           // if the genome is the same, order by the position in the genome
        return p1.second < p2.second;
    }
    
    return p1.first < p2.first;         // order by position of genomes in vector
}

bool compareGenomeMatch(GenomeMatch g1, GenomeMatch g2){
    if(g1.percentMatch == g2.percentMatch){
        return g1.genomeName < g2.genomeName;   // if percents are the same, order by name
    }
    
    return g1.percentMatch > g2.percentMatch;   // order by percents in descending order
}

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
    int m_minSearchLength;
    vector<Genome> m_genomes;
    
    // the Trie maps strings to pairs of ints       (position of genome in vector, position within genome)
    Trie<pair<int, int>> m_dna;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
:m_minSearchLength(minSearchLength) {}

int GenomeMatcherImpl::minimumSearchLength() const
{
    return m_minSearchLength;
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
    int pos = m_genomes.size();
    string frag = "";
    
    
    m_genomes.push_back(genome);
    
    for(int i = 0; i < genome.length() - m_minSearchLength + 1; i++){
        pair<int, int> p;
        p.first = pos;
        p.second = i;
        genome.extract(i, m_minSearchLength, frag);
        m_dna.insert(frag, p);
    }
}


bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    // return false for invalid input (lengths lower than minSearchLength)
    if(fragment.length() < minimumLength || minimumLength < m_minSearchLength){
        return false;
    }
    
    // get pairs with matching prefixes
    vector<pair<int, int>> dnaFragMatches;
    dnaFragMatches = m_dna.find(fragment.substr(0, m_minSearchLength), exactMatchOnly);
    
    // returns immdiately if there are no prefix matches
    if(dnaFragMatches.size() == 0)
        return false;
    
    // order the dna prefix matches by the genome they come from
    sort(dnaFragMatches.begin(), dnaFragMatches.end(), &(comparePairByGenome));
    
    // curGenome stores the position of the genome we are currently searching through
    int curGenome = dnaFragMatches[0].first;
    int longestPos = -1;
    int longest = -1;
    
    matches.clear();
    for(int i = 0; i < dnaFragMatches.size(); i++){
        // if curGenome and the next fragment's genome are different,
        // push the longest match (if there is a match) from the previous genome
        if(curGenome != dnaFragMatches[i].first){
            if(longest >= minimumLength){
                DNAMatch m;
                m.genomeName = m_genomes[curGenome].name();
                m.length = longest;
                m.position = longestPos;
                matches.push_back(m);
            }
            
            curGenome = dnaFragMatches[i].first;
            longest = -1;
            longestPos = -1;
        }
        
        bool errorsAllowed = !exactMatchOnly;
        int searchLength = fragment.length();
        string genomeFrag;
        
        // check whether extracting for the entire fragment's size would be out of bounds of the genome
        if(dnaFragMatches[i].second + fragment.length() > m_genomes[curGenome].length()){
            searchLength = m_genomes[curGenome].length() - dnaFragMatches[i].second;
        }
        
        // extract the most DNA bases as possible up to a max of fragment's length
        m_genomes[curGenome].extract(dnaFragMatches[i].second, searchLength, genomeFrag);

        // iterate until fragment and the extracted string aren't equal
        int curLength = 0;
        while(curLength < searchLength){
            if(fragment[curLength] != genomeFrag[curLength]){
                if(errorsAllowed){
                    errorsAllowed = false;
                    curLength++;
                    continue;
                }
                break;
            }
            curLength++;
        }
        
        // if the current fragment is a match and is longer than the previous match for this genome, update longestPos and longest
        if(curLength >= minimumLength && curLength > longest){
            longestPos = dnaFragMatches[i].second;
            longest = curLength;
        }
    }
    
    // add the longest match (if there is one) from the last genome
    if(longest >= minimumLength){
        DNAMatch m;
        m.genomeName = m_genomes[curGenome].name();
        m.length = longest;
        m.position = longestPos;
        matches.push_back(m);
    }
    
    if(matches.size() > 0)
        return true;
    return false;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    int numIterations = query.length()/fragmentMatchLength;
    map<string, int> numMatches;            // use a map to maintain the counts for the number of matches
    
    for(int i = 0; i < numIterations; i++){
        string frag;
        vector<DNAMatch> matches;
        
        query.extract(i*fragmentMatchLength, fragmentMatchLength, frag);
        findGenomesWithThisDNA(frag, fragmentMatchLength, exactMatchOnly, matches);
        
        for(int i = 0; i < matches.size(); i++){
            // first check if that genome already exists in the map
            map<string,int>::iterator it = numMatches.find(matches[i].genomeName);
            if(it == numMatches.end()){
                numMatches[matches[i].genomeName] = 1;
            }else{
                numMatches[matches[i].genomeName]++;
            }
        }
    }
    
    results.clear();
    for(map<string, int>::iterator it = numMatches.begin(); it != numMatches.end(); it++){
        double percent = (double)(*it).second/numIterations * 100 ; // percent is 0-100
        if(percent >= matchPercentThreshold){
            GenomeMatch g;
            g.genomeName = (*it).first;
            g.percentMatch = percent;
            
            results.push_back(g);
        }
    }
    sort(results.begin(), results.end(), compareGenomeMatch);
    
    if(results.size() > 0)
        return true;
    return false;
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
