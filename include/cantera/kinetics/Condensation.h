#ifndef CT_CONDENSATION_H
#define CT_CONDENSATION_H

#include "Kinetics.h"
#include "BulkKinetics.h"

namespace Cantera
{

class Condensation : public BulkKinetics {
 public:
  Condensation() = default;

  ~Condensation() override {}

  void setQuantityMoleFraction() {
    if (!m_use_mole_fraction) {
      m_ROP_ok = false;
    }
    m_use_mole_fraction = true;
  }

  void setQuantityConcentration() {
    if (m_use_mole_fraction) {
      m_ROP_ok = false;
    }
    m_use_mole_fraction = false;
  }

  bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;

  string kineticsType() const override {
    return "condensation";
  }

  bool isReversible(size_t i) override {
    return i < m_jxy.size() + m_jxxy.size() + m_jyy.size();
  }

  void updateROP() override;

 protected:
  //! reaction indices for x <=> y
  vector<size_t> m_jxy;

  //! reaction indices for x1 + x2 <=> y
  vector<size_t> m_jxxy;

  //! reaction indices for freezing reaction x(l) <=> x(s)
  vector<size_t> m_jyy;

  bool m_use_mole_fraction = false;
};

}

#endif
