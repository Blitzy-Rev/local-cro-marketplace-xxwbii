## WHY - Vision & Purpose

### Purpose & Users

We are solving the problem for small to mid-cap pharma users by providing an application that interfaces with an engine predicting chemical and molecular properties for small molecule drug discovery. The application will:

- Accept a CSV file containing a **SMILES column** and associated numerical property values.

- Allow users to **sort, filter, and organize** molecules into different libraries based on user-defined specifications.

- Enable users to **save, flag, and queue molecules** for different experimental workflows.

- Provide a **seamless workflow** between computational data aggregation and CRO (Contract Research Organization) services.

- Include an interface for CROs to receive molecule submissions, interact with users, set pricing, and return assay data.

The core users are **small pharma** managing molecular data and **CROs** providing wet lab testing services. This platform **streamlines the molecule-to-CRO pipeline**, making interactions efficient and organized.

----------

## WHAT - Core Requirements

### Functional Requirements

#### **User Management & Authentication**

- System must allow users to create and manage accounts with role-based access (**Pharma Users, CRO Users, Admins**).

- System must support secure login **without** external OAuth services or APIs.

#### **CSV Upload & Molecular Data Ingestion**

- Users must be able to upload a **CSV** containing:

  - A **SMILES column** (for molecular structures).

  - Multiple **numerical properties** (logP, MW, IC50, etc.).

- System must validate the CSV file format before ingestion.

- Users must be able to map CSV headers to **system-defined properties**.

- System must be adaptable to accept any new chemical property columns without requiring modifications.

#### **Molecule Sorting & Organization**

- Users must be able to **sort molecules** by property (e.g., MW, logP, affinity score).

- Users must be able to **filter molecules** based on numerical ranges, flags, and categories.

- Users must be able to **group molecules** into custom libraries (e.g., "High Binding Candidates", "Series 1").

- System must provide an **interactive drag-and-drop** interface for molecule organization.

#### **Molecule Management & Experiment Queuing**

- Users must be able to **flag molecules** for priority review.

- Users must be able to **add molecules to a queue** for experimental testing.

- Users must be able to **batch select molecules** and assign them to specific experimental runs.

- System must track molecule status (e.g., "Awaiting Testing", "Submitted to CRO", "Results Available").

#### **CRO Submission & Integration**

- Users must be able to:

  - **Select a CRO service** (e.g., binding assay, ADME profiling).

  - **Submit molecules** for testing via an interactive interface.

  - Attach **experimental specifications**, budget constraints, and legal documents.

- CRO Users must be able to:

  - Receive molecule submissions with attached metadata.

  - Communicate **assay pricing, turnaround times, and data requirements**.

  - Upload **assay results** and send them back to pharma users.

#### **User Interface & Experience (UI/UX)**

##### **Pharma User Interface**

- Intuitive **dashboard** for molecular data organization and experiment submission.

- **Drag-and-drop** functionality for molecule sorting and assignment.

- Clear visibility of **experiment status tracking** and CRO communications.

- Interactive tables with **real-time filtering and sorting** capabilities.

- **Auto-saving** of user configurations for faster workflow completion.

##### **CRO User Interface**

- **Submission review panel** for processing incoming molecular data.

- **Pricing and turnaround time estimator** for assays.

- **File upload functionality** for assay result submission.

- **Secure messaging system** for negotiations and clarifications with pharma users.

- **Legal and compliance document repository** for regulatory adherence.

##### **Overall UX Focus**

- **Minimalist and intuitive design** ensuring ease of use.

- **Seamless navigation** with clear step-by-step guidance.

- **Customizable views** for different user roles.

- **Real-time updates** for experiment progress and communication logs.

#### **Available Experimental Assessments**

The system must include a comprehensive list of experimental assessments, allowing users to select tests based on their research needs. Each test will request specific molecular data, such as SMILES strings and relevant properties. The platform will include the following categories:

##### **ADME (Absorption, Distribution, Metabolism, Excretion) Assays**

- Evaluate pharmacokinetic properties.

- Required data: SMILES, molecular weight, logP, solubility.

##### **Toxicity Assays**

- Assess cytotoxicity and mutagenicity potential.

- Required data: SMILES, hERG inhibition, cytotoxicity data.

##### **Enzyme Inhibition Assays**

- Determine inhibitory effects on target enzymes.

- Required data: SMILES, IC50, Ki values.

##### **Receptor Binding Assays**

- Measure binding affinity to receptors.

- Required data: SMILES, Kd values.

##### **Cell Viability Assays**

- Test compound effects on cell survival.

- Required data: SMILES, cytotoxicity data.

##### **Pharmacokinetics (PK) Studies**

- Examine absorption, distribution, metabolism, and excretion.

- Required data: SMILES, plasma concentration-time data.

##### **High-Throughput Screening (HTS) Assays**

- Screen large libraries for activity against targets.

- Required data: SMILES, fluorescence intensity.

##### **Genotoxicity Assays**

- Assess genetic mutation risks.

- Required data: SMILES, mutagenicity data.

##### **Metabolic Stability Assays**

- Evaluate compound breakdown by liver enzymes.

- Required data: SMILES, intrinsic clearance rates.

##### **Protein Binding Assays**

- Measure drug binding to plasma proteins.

- Required data: SMILES, protein binding percentages.

##### **CYP450 Inhibition Assays**

- Assess effects on cytochrome P450 enzymes.

- Required data: SMILES, inhibition constants (Ki).

##### **Transporter Interaction Assays**

- Evaluate drug interactions with transporters.

- Required data: SMILES, transporter substrate/inhibitor status.

##### **In Vivo Efficacy Studies**

- Test drug effects in animal models.

- Required data: SMILES, dosing regimen.

##### **Biomarker Analysis**

- Measure biological markers of drug effect.

- Required data: SMILES, biomarker data.

##### **Formulation Development**

- Optimize drug formulations.

- Required data: SMILES, solubility data.

##### **Stability Studies**

- Test compound stability under various conditions.

- Required data: SMILES, degradation profiles.

##### **Immunoassays**

- Detect and quantify proteins or antigens.

- Required data: SMILES, antigen/antibody info.

##### **Ion Channel Assays**

- Evaluate effects on ion channel function.

- Required data: SMILES, electrophysiological data.

##### **GPCR Assays**

- Assess interactions with GPCR drug targets.

- Required data: SMILES, receptor binding data.

These experimental assessments will be available in the platform’s marketplace, where users can request data, upload molecular information, and manage experiment submissions seamlessly.

----------

## HOW - Planning & Implementation

### **Deployment Strategy**

- System must be **fully deployable locally** using Docker containers.

- Deployment should be **streamlined** using a **multi-container Docker setup**.

- Local **Docker Desktop** should be sufficient for setting up the full application.

- No cloud or external hosting dependencies should be required for full functionality.

----------

## **Final Notes**

This system must be **completely independent**, relying on **no external APIs** for any core functions. All functionality, from authentication to AI processing, will be handled **locally**. Deployment will be simplified using **Docker**, ensuring seamless installation and operation on any local machine.