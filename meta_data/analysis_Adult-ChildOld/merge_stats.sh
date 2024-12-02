for file in `ls multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '6p' "$file" >> mean_stats.txt
done

for file in `ls multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '8p' "$file" >> mean_stats.txt
done

for file in `ls multivoxel_extract_mean_good_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '9p' "$file" >> mean_stats.txt
done



for file in `ls multivoxel_extract_MyLinearModel_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '6p' "$file" >> contrast_stats.txt
done

for file in `ls multivoxel_extract_MyLinearModel_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '8p' "$file" >> contrast_stats.txt
done

for file in `ls multivoxel_extract_MyLinearModel_z_voxelCorrected_p_0.00100_10_blob*.txt` ; do
    sed -n '9p' "$file" >> contrast_stats.txt
done