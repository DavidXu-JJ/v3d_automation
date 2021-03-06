function save_v3d_pointcloud_file(m_struct, filename)
%function save_v3d_pointcloud_file(m_struct, filename)
%
% Save the .apo point cloud format data file used in V3D
% 
% m_struct will consist of point cloud coordinates and other information
% (e.g. name/comments/types/shape)
%
% V3D website: see software page of http://penglab.janelia.org
%
% by Hanchuan Peng
% 20090724

fp = fopen(filename, 'w');
if (fp<0),
    disp('Fail to open the file to save V3D pointcloud (.apo) format data');
    return;
end;

fprintf(fp, ...
    ['#no, order_info, name, comment, z, x, y, ' ...
    'max_intensity, mean_intensity, sdev_intensity, ' ...
    'rgn size (#voxels), mass, reserved_anno_1, reserved_anno_2, ' ...
    'reserved_anno_3, color_red, color_green, color_blue\n']);

for i=1:length(m_struct),
  S = m_struct{i};  
  fprintf(fp, '%d, %s, %s, %s, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %s, %s, %s, %d, %d, %d\n', ...
      S.n, ...
      trimmed_str(S.orderinfo, 99), ...
      trimmed_str(S.name, 99), ...
      trimmed_str(S.comment, 99), ...
      S.z, ...
      S.x, ...
      S.y, ...
      S.pixmax, ...
      S.intensity, ...
      S.sdev, ...
      S.volsize, ...
      S.mass, ...
      '', ...
      '', ...
      '', ...
      S.color.r, ...
      S.color.g, ...
      S.color.b); 
end;

fclose(fp);

return;


