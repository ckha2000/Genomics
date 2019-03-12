#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
    string m_name;
    string m_genome;
    int m_length;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
:m_name(nm), m_genome(sequence), m_length(sequence.length())
{}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
    if(!genomeSource)
        return false;
    
    string temp, tempName, tempGenome;
    getline(genomeSource, temp);
    
    if(temp[0] != '>')
        return false;
    
    tempName = temp.substr(1);      // get the name of the genome
    tempGenome = "";
    
    while(genomeSource){
        getline(genomeSource, temp);
        
        for(int i = 0; i < temp.length(); i++){
            if(temp[i] == '>'){                 // add a new genome every time we reach a name line
                Genome g = Genome(tempName, tempGenome);
                genomes.push_back(g);
                if(temp.length() == 1)          // improper format if '>' not followed by a valid name
                    return false;

                tempName = temp.substr(1);
                tempGenome = "";
                break;
            }
            switch(temp[i]){
                case 'g':
                case 'G':
                    tempGenome += 'G';
                    break;
                case 'a':
                case 'A':
                    tempGenome += 'A';
                    break;
                case 'c':
                case 'C':
                    tempGenome += 'C';
                    break;
                case 't':
                case 'T':
                    tempGenome += 'T';
                    break;
                case 'n':
                case 'N':
                    tempGenome += 'N';
                    break;
                default:
                    return false;
            }
        }
    }
    if(tempGenome != ""){          // if there are genome lines after the last name line, add it to the vector
        Genome g = Genome(tempName, tempGenome);
        genomes.push_back(g);
        return true;
    }
    return false;                   // if there are no genome lines, return false
}

int GenomeImpl::length() const
{
    return m_length;
}

string GenomeImpl::name() const
{
    return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
    if(position + length > m_length){
        return false;
    }
    fragment = m_genome.substr(position, length);
    return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
