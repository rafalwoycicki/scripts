python
import sys
import numpy as np
import doubletdetection
import numba.targets
import tarfile
import matplotlib.pyplot as plt
##################################################

Nod1_matrix_path = '../results/Nod1_filtered_feature_bc_matrix/matrix.mtx.gz'
Nod1_raw_counts = doubletdetection.load_mtx(Nod1_matrix_path)
Nod1_zero_genes = (np.sum(Nod1_raw_counts, axis=0) == 0).A.ravel()
Nod1_raw_counts = Nod1_raw_counts[:, ~Nod1_zero_genes]
Nod1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
Nod1_old_doublets = Nod1_clf_old.fit(Nod1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

Nod1_db_old = Nod1_old_doublets.nonzero()
f = open("./Nod1_ddold_prepathms.txt", "w")
for d in Nod1_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
###########################
Nod2_matrix_path = '../results/Nod2_filtered_feature_bc_matrix/matrix.mtx.gz'
Nod2_raw_counts = doubletdetection.load_mtx(Nod2_matrix_path)
Nod2_zero_genes = (np.sum(Nod2_raw_counts, axis=0) == 0).A.ravel()
Nod2_raw_counts = Nod2_raw_counts[:, ~Nod2_zero_genes]
Nod2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
Nod2_old_doublets = Nod2_clf_old.fit(Nod2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

Nod2_db_old = Nod2_old_doublets.nonzero()

f = open("./Nod2_ddold_prepathms.txt", "w")
for d in Nod2_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
################
#############################Make NodsBdiaz
Nod1Bdiaz_matrix_path = '../results/Nod1_v4nuc_BdiazUSDA/filtered_feature_bc_matrix/matrix.mtx.gz'
Nod1Bdiaz_raw_counts = doubletdetection.load_mtx(Nod1Bdiaz_matrix_path)
Nod1Bdiaz_zero_genes = (np.sum(Nod1Bdiaz_raw_counts, axis=0) == 0).A.ravel()
Nod1Bdiaz_raw_counts = Nod1Bdiaz_raw_counts[:, ~Nod1Bdiaz_zero_genes]
Nod1Bdiaz_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
Nod1Bdiaz_old_doublets = Nod1Bdiaz_clf_old.fit(Nod1Bdiaz_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

Nod1Bdiaz_db_old = Nod1Bdiaz_old_doublets.nonzero()
f = open("./Nod1Bdiaz_ddold_prepathms.txt", "w")
for d in Nod1Bdiaz_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
###########################
Nod2Bdiaz_matrix_path = '../results/Nod2_v4nuc_BdiazUSDA/filtered_feature_bc_matrix/matrix.mtx.gz'
Nod2Bdiaz_raw_counts = doubletdetection.load_mtx(Nod2Bdiaz_matrix_path)
Nod2Bdiaz_zero_genes = (np.sum(Nod2Bdiaz_raw_counts, axis=0) == 0).A.ravel()
Nod2Bdiaz_raw_counts = Nod2Bdiaz_raw_counts[:, ~Nod2Bdiaz_zero_genes]
Nod2Bdiaz_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
Nod2Bdiaz_old_doublets = Nod2Bdiaz_clf_old.fit(Nod2Bdiaz_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

Nod2Bdiaz_db_old = Nod2Bdiaz_old_doublets.nonzero()

f = open("./Nod2Bdiaz_ddold_prepathms.txt", "w")
for d in Nod2Bdiaz_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()



###############
TL1_matrix_path = '../results/TL1_filtered_feature_bc_matrix/matrix.mtx.gz'
TL1_raw_counts = doubletdetection.load_mtx(TL1_matrix_path)
TL1_zero_genes = (np.sum(TL1_raw_counts, axis=0) == 0).A.ravel()
TL1_raw_counts = TL1_raw_counts[:, ~TL1_zero_genes]
TL1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
TL1_old_doublets = TL1_clf_old.fit(TL1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

TL1_db_old = TL1_old_doublets.nonzero()

f = open("./TL1_ddold_prepathms.txt", "w")
for d in TL1_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
#################
Tri_matrix_path = '../results/Tri_filtered_feature_bc_matrix/matrix.mtx.gz'
Tri_raw_counts = doubletdetection.load_mtx(Tri_matrix_path)
Tri_zero_genes = (np.sum(Tri_raw_counts, axis=0) == 0).A.ravel()
Tri_raw_counts = Tri_raw_counts[:, ~Tri_zero_genes]
Tri_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
Tri_old_doublets = Tri_clf_old.fit(Tri_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

Tri_db_old = Tri_old_doublets.nonzero()

f = open("./Tri_ddold_prepathms.txt", "w")
for d in Tri_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
#####################################################################################

tri2_matrix_path = '../results/tri2_filtered_feature_bc_matrix/matrix.mtx.gz'
tri2_raw_counts = doubletdetection.load_mtx(tri2_matrix_path)
tri2_zero_genes = (np.sum(tri2_raw_counts, axis=0) == 0).A.ravel()
tri2_raw_counts = tri2_raw_counts[:, ~tri2_zero_genes]
tri2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
tri2_old_doublets = tri2_clf_old.fit(tri2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

tri2_db_old = tri2_old_doublets.nonzero()

f = open("./tri2_ddold_prepathms.txt", "w")
for d in tri2_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
##################
############################################################################################
########################################################################################
RootMock_matrix_path = '../results/RootMock_filtered_feature_bc_matrix/matrix.mtx.gz'
RootMock_raw_counts = doubletdetection.load_mtx(RootMock_matrix_path)
RootMock_zero_genes = (np.sum(RootMock_raw_counts, axis=0) == 0).A.ravel()
RootMock_raw_counts = RootMock_raw_counts[:, ~RootMock_zero_genes]
RootMock_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
RootMock_old_doublets = RootMock_clf_old.fit(RootMock_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

RootMock_db_old = RootMock_old_doublets.nonzero()

f = open("./RootMock_ddold_prepathms.txt", "w")
for d in RootMock_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
##############
RootBjap_matrix_path = '../results/RootBjap_filtered_feature_bc_matrix/matrix.mtx.gz'
RootBjap_raw_counts = doubletdetection.load_mtx(RootBjap_matrix_path)
RootBjap_zero_genes = (np.sum(RootBjap_raw_counts, axis=0) == 0).A.ravel()
RootBjap_raw_counts = RootBjap_raw_counts[:, ~RootBjap_zero_genes]
RootBjap_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
RootBjap_old_doublets = RootBjap_clf_old.fit(RootBjap_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

RootBjap_db_old = RootBjap_old_doublets.nonzero()

f = open("./RootBjap_ddold_prepathms.txt", "w")
for d in RootBjap_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
##############
LM2_matrix_path = '../results/LM2_filtered_feature_bc_matrix/matrix.mtx.gz'
LM2_raw_counts = doubletdetection.load_mtx(LM2_matrix_path)
LM2_zero_genes = (np.sum(LM2_raw_counts, axis=0) == 0).A.ravel()
LM2_raw_counts = LM2_raw_counts[:, ~LM2_zero_genes]
LM2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
LM2_old_doublets = LM2_clf_old.fit(LM2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

LM2_db_old = LM2_old_doublets.nonzero()

f = open("./LM2_ddold_prepathms.txt", "w")
for d in LM2_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
###############
LM1_matrix_path = '../results/LM1_filtered_feature_bc_matrix/matrix.mtx.gz'
LM1_raw_counts = doubletdetection.load_mtx(LM1_matrix_path)
LM1_zero_genes = (np.sum(LM1_raw_counts, axis=0) == 0).A.ravel()
LM1_raw_counts = LM1_raw_counts[:, ~LM1_zero_genes]
LM1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
LM1_old_doublets = LM1_clf_old.fit(LM1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

LM1_db_old = LM1_old_doublets.nonzero()

f = open("./LM1_ddold_prepathms.txt", "w")
for d in LM1_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
################
FLOW_matrix_path = '../results/FLOW_filtered_feature_bc_matrix/matrix.mtx.gz'
FLOW_raw_counts = doubletdetection.load_mtx(FLOW_matrix_path)
FLOW_zero_genes = (np.sum(FLOW_raw_counts, axis=0) == 0).A.ravel()
FLOW_raw_counts = FLOW_raw_counts[:, ~FLOW_zero_genes]
FLOW_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
FLOW_old_doublets = FLOW_clf_old.fit(FLOW_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

FLOW_db_old = FLOW_old_doublets.nonzero()

f = open("./FLOW_ddold_prepathms.txt", "w")
for d in FLOW_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
##############
FLOW2_matrix_path = '../results/FLOW2_filtered_feature_bc_matrix/matrix.mtx.gz'
FLOW2_raw_counts = doubletdetection.load_mtx(FLOW2_matrix_path)
FLOW2_zero_genes = (np.sum(FLOW2_raw_counts, axis=0) == 0).A.ravel()
FLOW2_raw_counts = FLOW2_raw_counts[:, ~FLOW2_zero_genes]
FLOW2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
FLOW2_old_doublets = FLOW2_clf_old.fit(FLOW2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

FLOW2_db_old = FLOW2_old_doublets.nonzero()

f = open("./FLOW2_ddold_prepathms.txt", "w")
for d in FLOW2_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
#######

EM2_matrix_path = '../results/EM2_filtered_feature_bc_matrix/matrix.mtx.gz'
EM2_raw_counts = doubletdetection.load_mtx(EM2_matrix_path)
EM2_zero_genes = (np.sum(EM2_raw_counts, axis=0) == 0).A.ravel()
EM2_raw_counts = EM2_raw_counts[:, ~EM2_zero_genes]
EM2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
EM2_old_doublets = EM2_clf_old.fit(EM2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

EM2_db_old = EM2_old_doublets.nonzero()

f = open("./EM2_ddold_prepathms.txt", "w")
for d in EM2_db_old[0] :
 f.write(str(d+1) + "\n")
#[enter] - here you physically press ENTER/RETURN key
f.close()

##############
EM1_matrix_path = '../results/EM1_filtered_feature_bc_matrix/matrix.mtx.gz'
EM1_raw_counts = doubletdetection.load_mtx(EM1_matrix_path)
EM1_zero_genes = (np.sum(EM1_raw_counts, axis=0) == 0).A.ravel()
EM1_raw_counts = EM1_raw_counts[:, ~EM1_zero_genes]
EM1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
EM1_old_doublets = EM1_clf_old.fit(EM1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

EM1_db_old = EM1_old_doublets.nonzero()

f = open("./EM1_ddold_prepathms.txt", "w")
for d in EM1_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()
####################





################################################

B1_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyBjap1_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path
B2_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyBjap2_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path
B3_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyBjap3_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path
M1_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyMock1_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path
M2_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyMock2_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path
M3_matrix_path = '/home/ravwoy/WORK2/LibaultM/Project1_Soy/SoyMock3_v4_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz' # change path

B1_raw_counts = doubletdetection.load_mtx(B1_matrix_path) 
B2_raw_counts = doubletdetection.load_mtx(B2_matrix_path)
B3_raw_counts = doubletdetection.load_mtx(B3_matrix_path) 
M1_raw_counts = doubletdetection.load_mtx(M1_matrix_path)
M2_raw_counts = doubletdetection.load_mtx(M2_matrix_path) 
M3_raw_counts = doubletdetection.load_mtx(M3_matrix_path)

B1_zero_genes = (np.sum(B1_raw_counts, axis=0) == 0).A.ravel()
B2_zero_genes = (np.sum(B2_raw_counts, axis=0) == 0).A.ravel()
B3_zero_genes = (np.sum(B3_raw_counts, axis=0) == 0).A.ravel()
M1_zero_genes = (np.sum(M1_raw_counts, axis=0) == 0).A.ravel()
M2_zero_genes = (np.sum(M2_raw_counts, axis=0) == 0).A.ravel()
M3_zero_genes = (np.sum(M3_raw_counts, axis=0) == 0).A.ravel()

B1_raw_counts = B1_raw_counts[:, ~B1_zero_genes]
B2_raw_counts = B2_raw_counts[:, ~B2_zero_genes]
B3_raw_counts = B3_raw_counts[:, ~B3_zero_genes]
M1_raw_counts = M1_raw_counts[:, ~M1_zero_genes]
M2_raw_counts = M2_raw_counts[:, ~M2_zero_genes]
M3_raw_counts = M3_raw_counts[:, ~M3_zero_genes]

#clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
#L5_doublets = clf.fit(L5_raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)
#L6_doublets = clf.fit(L6_raw_counts).predict(p_thresh=1e-16, voter_thresh=0.5)

#################################################

B1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
B1_old_doublets = B1_clf_old.fit(B1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

B1_db_old = B1_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/B1_ddold_prepathms.txt", "w")
for d in B1_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()

################################################

B2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
B2_old_doublets = clf_old.fit(B2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

B2_db_old = B2_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/B2_ddold_prepathms.txt", "w")
for d in B2_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()

#######################################################

B3_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
B3_old_doublets = clf_old.fit(B3_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

B3_db_old = B3_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/B3_ddold_prepathms.txt", "w")
for d in B3_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()

##############################################################

M1_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
M1_old_doublets = clf_old.fit(M1_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

M1_db_old = M1_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/M1_ddold_prepathms.txt", "w")
for d in M1_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()

#################################################################

M2_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
M2_old_doublets = clf_old.fit(M2_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

M2_db_old = M2_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/M2_ddold_prepathms.txt", "w")
for d in M2_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()

########################################################################

M3_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
M3_old_doublets = clf_old.fit(M3_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

M3_db_old = M3_old_doublets.nonzero()

f = open("/home/ravwoy/WORK2/LibaultM/Project1_Soy/M3_ddold_prepathms.txt", "w")
for d in M3_db_old[0] :
 f.write(str(d+1) + "\n")

f.close()
#############

#L5_db = L5_doublets.nonzero()
#L6_db = L6_doublets.nonzero()

###

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SoyMock_Medicago_matrix_path = './soy4_and_medicago_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
SoyMock_Medicago_raw_counts = doubletdetection.load_mtx(SoyMock_Medicago_matrix_path)
SoyMock_Medicago_zero_genes = (np.sum(SoyMock_Medicago_raw_counts, axis=0) == 0).A.ravel()
SoyMock_Medicago_raw_counts = SoyMock_Medicago_raw_counts[:, ~SoyMock_Medicago_zero_genes]
SoyMock_Medicago_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
SoyMock_Medicago_old_doublets = SoyMock_Medicago_clf_old.fit(SoyMock_Medicago_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

SoyMock_Medicago_db_old = SoyMock_Medicago_old_doublets.nonzero()
f = open("./SoyMock_Medicago_ddold_prepathms.txt", "w")
for d in SoyMock_Medicago_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SoyBjap_Corn_matrix_path = './soy4_and_may5_count/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
SoyBjap_Corn_raw_counts = doubletdetection.load_mtx(SoyBjap_Corn_matrix_path)
SoyBjap_Corn_zero_genes = (np.sum(SoyBjap_Corn_raw_counts, axis=0) == 0).A.ravel()
SoyBjap_Corn_raw_counts = SoyBjap_Corn_raw_counts[:, ~SoyBjap_Corn_zero_genes]
SoyBjap_Corn_clf_old = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=True, standard_scaling=False, verbose=True)
SoyBjap_Corn_old_doublets = SoyBjap_Corn_clf_old.fit(SoyBjap_Corn_raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

SoyBjap_Corn_db_old = SoyBjap_Corn_old_doublets.nonzero()
f = open("./SoyBjap_Corn_ddold_prepathms.txt", "w")
for d in SoyBjap_Corn_db_old[0] :
 f.write(str(d+1) + "\n")
[enter]
f.close()

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
