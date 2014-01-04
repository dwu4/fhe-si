#include "Plaintext.h"
#include "NumbTh.h"
#include <algorithm>
#include "assert.h"

void Plaintext::Init(const ZZ_pX &msg) {
  message = msg;
}

void Plaintext::Init(const vector<ZZ_pX> &msgs) {
  EmbedInSlots(msgs);
}

Plaintext &Plaintext::operator=(const Plaintext &other) {
  assert(&context == &other.context);
  message = other.message;
  
  return *this;
}

bool Plaintext::operator==(const Plaintext &other) const {
  return ((&context == &other.context) && 
          (message == other.message));
}

void Plaintext::EmbedInSlots(const vector<ZZ_pX> &msgs, bool onlyUsable) {
  context.GetPlaintextSpace().EmbedInSlots(message, msgs, onlyUsable);
}

void Plaintext::DecodeSlots(vector<ZZ_pX> &msgBatch, bool onlyUsable) {
  context.GetPlaintextSpace().DecodeSlots(msgBatch, message, onlyUsable);
}

ostream &operator<<(ostream &os, const Plaintext &ptxt) {
  return (os << ptxt.message);
}