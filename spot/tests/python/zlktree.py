# -*- mode: python; coding: utf-8 -*-
# Copyright (C) 2021 Laboratoire de Recherche et Développement de l'Epita
# (LRDE).
#
# This file is part of Spot, a model checking library.
#
# Spot is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Spot is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import spot

a = spot.automaton("""HOA: v1 States: 5 Start: 0 AP: 2 "p0" "p1"
Acceptance: 4 Inf(3) | Fin(3) properties: trans-labels explicit-labels
trans-acc --BODY-- State: 0 [!0&!1] 3 [!0&!1] 4 State: 1 [!0&!1] 4 {3}
[0&!1] 0 {2} [!0&1] 1 {2} State: 2 [!0&1] 0 {0 2} [!0&!1] 1 State: 3
[!0&1] 2 State: 4 [0&!1] 3 --END--""")
b = spot.zielonka_tree_transform(a)
assert spot.are_equivalent(a, b)
assert b.acc().is_buchi()

def report_missing_exception():
    raise RuntimeError("missing exception")

a = spot.automaton("""
HOA: v1 States: 10 Start: 0 AP: 2 "p0" "p1" acc-name: Rabin 3
Acceptance: 6 (Fin(0) & Inf(1)) | (Fin(2) & Inf(3)) | (Fin(4) &
Inf(5)) properties: trans-labels explicit-labels trans-acc --BODY--
State: 0 [0&!1] 0 {1 2 3} [!0&!1] 8 [0&1] 4 State: 1 [0&1] 1 {2}
State: 2 [0&1] 8 {3} [0&1] 2 {1} [!0&1] 4 {3 4} [!0&!1] 3 {2 5} State:
3 [!0&!1] 5 [0&1] 8 {1 2} [!0&!1] 9 {1} State: 4 [0&!1] 3 {0 2} [!0&1]
5 {1} State: 5 [!0&!1] 6 [!0&!1] 7 {2} [0&!1] 3 {1} [!0&1] 5 State: 6
[!0&!1] 1 [!0&!1] 2 {4} [0&!1] 0 {1 3 4} [0&1] 5 [!0&!1] 3 State: 7
[0&1] 4 {0} [0&1] 5 [0&1] 0 {3} [0&1] 1 {2 4} State: 8 [0&!1] 6 {0 4
5} [!0&!1] 7 {4} [0&!1] 2 {1 3} [0&1] 0 {0 1 4} State: 9 [!0&1] 6 {4}
[!0&!1] 2 {5} [!0&!1] 0 {3} [!0&!1] 5 --END--""")
aa = spot.acd(a)
try:
    assert aa.has_rabin_shape()
except RuntimeError as e:
    assert 'CHECK_RABIN' in str(e)
else:
    report_missing_exception()

try:
    assert not aa.has_streett_shape()
except RuntimeError as e:
    assert 'CHECK_STREETT' in str(e)
else:
    report_missing_exception()

try:
    assert not aa.has_parity_shape()
except RuntimeError as e:
    assert 'CHECK_PARITY' in str(e)
else:
    report_missing_exception()


aa = spot.acd(a, spot.acd_options_CHECK_RABIN)
assert aa.has_rabin_shape()
assert aa.node_count() == 13

try:
    assert not aa.has_streett_shape()
except RuntimeError as e:
    assert 'CHECK_STREETT' in str(e)
else:
    report_missing_exception()

try:
    assert aa.has_parity_shape()
except RuntimeError as e:
    assert 'CHECK_PARITY' in str(e)
else:
    report_missing_exception()

aa = spot.acd(a, (spot.acd_options_CHECK_PARITY |
                  spot.acd_options_ABORT_WRONG_SHAPE))
assert aa.has_rabin_shape()
assert not aa.has_streett_shape()
assert not aa.has_parity_shape()
assert aa.node_count() == 0
try:
    aa.first_branch(0)
except RuntimeError as e:
    assert 'ABORT_WRONG_SHAPE' in str(e)
else:
    report_missing_exception()

try:
    aa.step(0, 0)
except RuntimeError as e:
    assert 'incorrect branch number' in str(e)
else:
    report_missing_exception()

try:
    aa.node_acceptance(0)
except RuntimeError as e:
    assert 'unknown node' in str(e)
else:
    report_missing_exception()

try:
    aa.edges_of_node(0)
except RuntimeError as e:
    assert 'unknown node' in str(e)
else:
    report_missing_exception()

try:
    aa.node_level(0)
except RuntimeError as e:
    assert 'unknown node' in str(e)
else:
    report_missing_exception()

a = spot.translate('true')
a.set_acceptance(spot.acc_cond('f'))
b = spot.acd_transform(a)
assert a.equivalent_to(b)

a = spot.translate('true')
a.set_acceptance(spot.acc_cond('f'))
b = spot.zielonka_tree_transform(a)
assert a.equivalent_to(b)

a = spot.automaton("""HOA: v1 name: "^ G F p0 G F p1" States: 5 Start:
2 AP: 2 "a" "b" acc-name: Rabin 2 Acceptance: 4 (Fin(0) & Inf(1)) |
(Fin(2) & Inf(3)) properties: trans-labels explicit-labels state-acc
complete properties: deterministic --BODY-- State: 0 {0} [!0&!1] 0
[0&!1] 4 [!0&1] 3 [0&1] 2 State: 1 {2} [!0&!1] 1 [0&!1] 4 [!0&1] 3
[0&1] 2 State: 2 {0 2} [!0&!1] 2 [0&!1] 4 [!0&1] 3 [0&1] 2 State: 3 {1
2} [!0&!1] 1 [0&!1] 4 [!0&1] 3 [0&1] 2 State: 4 {0 3} [!0&!1] 0 [0&!1]
4 [!0&1] 3 [0&1] 2 --END--""")
b = spot.acd_transform_sbacc(a, True)
assert str(b.acc()) == '(3, Fin(0) & (Inf(1) | Fin(2)))'
assert a.equivalent_to(b)
b = spot.acd_transform_sbacc(a, False)
assert str(b.acc()) == '(2, Fin(0) & Inf(1))'
assert a.equivalent_to(b)
