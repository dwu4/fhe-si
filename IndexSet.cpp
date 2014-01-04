/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
 
#include "IndexSet.h"

const IndexSet& IndexSet::emptySet()
{
   static IndexSet empty;
   return empty;
}

  // constructs an interval, low to high
void IndexSet::intervalConstructor(long low, long high) {
  assert(low >= 0);

  if (high < low) {
    _first = 0; _last = -1; _card = 0;
  }
  else {
    rep.resize(high+1);
    for (long i = 0; i < low; i++)
      rep[i] = false;
    for (long i = low; i <= high; i++)
      rep[i] = true;

    _first = low; _last = high; _card = high-low+1;
  }
}

long IndexSet::next(long j) const {
  if (_card == 0) return j+1;
  if (j >= _last) return j + 1;
  if (j < _first) return _first;
  j++;
  while (rep[j] == false) j++;
  return j;
}

long IndexSet::prev(long j) const {
  if (_card == 0) return j-1;
  if (j > _last) return _last;
  if (j <= _first) return j-1;
  j--;
  while (rep[j] == false) j--;
  return j;
}

bool IndexSet::contains(long j) const {
  if (j < _first || j > _last) return false;
  return rep[j];
}


bool IndexSet::contains(const IndexSet& s) const {
  for (long i = s.first(); i <= s.last(); i = s.next(i))
    if (!contains(i)) return false;
  return true;
}


bool IndexSet::disjointFrom(const IndexSet& s) const
{
  // quick tests for some common cases
  if (card() == 0 || s.card() == 0
      || last() < s.first() || s.last() < first()) return true;

  for (long i = s.first(); i <= s.last(); i = s.next(i))
    if (contains(i)) return false;
  return true;
}



bool IndexSet::operator==(const IndexSet& s) const {
  if (this == &s) return true;
  if (_card != s._card) return false;
  if (_first != s._first) return false;
  if (_last != s._last) return false;

  return equal(rep.begin()+_first, rep.begin()+_last+1, 
               s.rep.begin()+_first); 
  // NOTE: maybe vector<bool> optimizes this???
}

void IndexSet::clear() {
  rep.resize(0);
  _first = 0; _last = -1; _card = 0;
}

void IndexSet::insert(long j) {
  assert(j >= 0);

  long oldSize = rep.size();
  if (j >= oldSize) {
    rep.resize(j+1);
    for (long i = oldSize; i <= j; i++) rep[i] = false;
  }

  if (_card == 0) {
    _first = _last = j;
    _card = 1;
  }
  else {
    if (j > _last) _last = j;
    if (j < _first) _first = j;
    if (rep[j] == false) _card++;
  }

  rep[j] = true;
}

void IndexSet::remove(long j) {
  assert(j >= 0);

  if (j >= (long) rep.size()) return;
  if (rep[j] == false) return;

  long newFirst = _first, newLast = _last;

  if (_card == 1) {
    newFirst = 0;
    newLast = -1;
  }
  else {
    if (_last == j) newLast = prev(_last);
    if (_first == j) newFirst = next(_first);
  }

  _first = newFirst;
  _last = newLast;
  _card--;
  rep[j] = false;
}

void IndexSet::insert(const IndexSet& s) {
  if (this == &s) return;
  if (s.card() == 0) return;
  if (card() == 0) {
    *this = s;
    return;
  }

  for (long i = s.last(); i >= s.first(); i = s.prev(i)) insert(i);
  // NOTE: traversal done from high to low so as to trigger at 
  // at most one resize

}

void IndexSet::remove(const IndexSet& s) {
  if (this == &s) { clear(); return; }
  if (s.card() == 0) return;
  if (card() == 0) return;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) remove(i);
  // NOTE: traversal order should not matter here
}


void IndexSet::retain(const IndexSet& s) {
  if (this == &s) return;
  if (s.card() == 0) { clear(); return; }
  if (card() == 0) return;

  for (long i = first(); i <= last(); i = next(i)) {
    if (!s.contains(i)) remove(i);
  }
}

// union
IndexSet operator|(const IndexSet& s, const IndexSet& t) {
  IndexSet r = s;
  r.insert(t);
  return r;
}

// intersection
IndexSet operator&(const IndexSet& s, const IndexSet& t) {
  IndexSet r = s;
  r.retain(t);
  return r;
}

// exclusive-or
IndexSet operator^(const IndexSet& s, const IndexSet& t) {
  IndexSet r = s | t;
  r.remove(s & t);
  return r;
}

// set minus
IndexSet operator/(const IndexSet& s, const IndexSet& t) {
  IndexSet r = s;
  r.remove(t);
  return r;
}


ostream& operator << (ostream& str, const IndexSet& set)
{
  if (set.card() == 0) {
    str << "[]";
  }
  else if (set.card() == 1) {
    str << "[" << set.first() << "]";
  }
  else {
    str << "[" << set.first();
    for (long i = set.next(set.first()); i <= set.last(); i = set.next(i))
      str << " " << i;
    str << "]";
  }

  return str;
}


// functional card
long card(const IndexSet& s) { return s.card(); }

// functional "contains"
bool operator<=(const IndexSet& s1, const IndexSet& s2) {
  return s2.contains(s1);
}

bool operator<(const IndexSet& s1, const IndexSet& s2) {
  return card(s1) < card(s2) && s2.contains(s1);
}

bool operator>=(const IndexSet& s1, const IndexSet& s2) {
  return s1.contains(s2);
}

bool operator>(const IndexSet& s1, const IndexSet& s2) {
  return card(s2) < card(s1) && s1.contains(s2);
}


#if 0
int main()
{
  cout << "hello world\n";

  IndexSet s1 = IndexSet(1) | IndexSet(3) | IndexSet(5);

  IndexSet s2 = IndexSet(2, 5);

  s1 = s1 / s2;
 
  cout << s1.first() << "\n";
  cout << s1.last() << "\n";
  cout << s1.card() << "\n";

  cout << s1 << "\n";
}
#endif
       

