import pydicom
import os
import csv

# ----------------------------
# CONFIGURATION
# ----------------------------
root_dir = '/storage/group/alh98/default/VLN_BIDS/participant_data_complete'  # top-level folder with subjects
output_csv = '/storage/work/krf5429/GNG-Inhibition-Task/Data/participants.csv'  # full path for CSV
max_subjects = 34  # only process first 34 folders
task_keyword = 'GO_NOGO_RUN1'  # folder name pattern for the task run

# ----------------------------
# SCRIPT
# ----------------------------
subjects = {}

# Get first 34 subject folders
all_subj_folders = sorted([f for f in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, f))])[:max_subjects]

for subj_folder in all_subj_folders:
    subj_path = os.path.join(root_dir, subj_folder)
    session1_path = os.path.join(subj_path, 'session1')

    if not os.path.isdir(session1_path):
        print(f"No session1 folder for subject {subj_folder}, skipping")
        continue

    # Find the GO/No-Go RUN1 folder
    task_folder = None
    for f in os.listdir(session1_path):
        if os.path.isdir(os.path.join(session1_path, f)) and task_keyword in f.upper():
            task_folder = f
            break

    if task_folder is None:
        print(f"No GO/No-Go RUN1 folder found for subject {subj_folder}, skipping")
        continue

    task_path = os.path.join(session1_path, task_folder)

    # Get first DICOM file
    dicom_files = sorted([f for f in os.listdir(task_path) if f.lower().endswith('.dcm')])
    if not dicom_files:
        print(f"No DICOM files found in {task_path}, skipping")
        continue

    first_dcm = os.path.join(task_path, dicom_files[0])

    # Read DICOM metadata
    try:
        ds = pydicom.dcmread(first_dcm, stop_before_pixels=True)
        dicom_pid = getattr(ds, "PatientID", "MISSING")
        sex = getattr(ds, "PatientSex", "MISSING")
        birthdate = getattr(ds, "PatientBirthDate", "MISSING")
        age_str = getattr(ds, "PatientAge", "MISSING")
        # Convert age '030Y' -> 30
        age_num = int(age_str[:-1]) if age_str != "MISSING" and age_str.endswith('Y') else "MISSING"

        subjects[subj_folder] = {
            "SubjectID": subj_folder,   # folder name
            "PatientID": dicom_pid,     # DICOM PatientID
            "Sex": sex,
            "Age": age_num,
            "BirthDate": birthdate,
            "TaskFolder": task_folder
        }

        print(f"Processed subject {subj_folder}: PatientID={dicom_pid}, Age={age_num}, Sex={sex}")

    except Exception as e:
        print(f"Failed to read DICOM for subject {subj_folder}: {e}")

# ----------------------------
# WRITE CSV
# ----------------------------
# Make sure the output directory exists
os.makedirs(os.path.dirname(output_csv), exist_ok=True)

with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ["SubjectID", "PatientID", "Sex", "Age", "BirthDate", "TaskFolder"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for sub in subjects.values():
        writer.writerow(sub)

print(f"CSV file created: {output_csv}")