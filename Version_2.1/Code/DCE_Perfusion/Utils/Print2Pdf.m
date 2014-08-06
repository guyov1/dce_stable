function [fig_index_updated] = Print2Pdf(fig_handle, fig_index, fig_name, out_directory, title, subtitle)
% Print2Pdf - printing figure to pdf 

            fig_file = [out_directory fig_name];
            gprint(fig_handle,fig_file);
            %gprint(fig_handle,'Run_Output/Larsson_vs_Sourbron.png');
            
            % If title is needed
            if (~strcmp(title,''))
                fig_index    = fig_index + 1;
                idx_string   = ['idx_' num2str(fig_index,'%03i')];
                %AddToLog('Run_Output/',idx_string,'\\subsection*{\\underline{Larsson Vs. Sourbron filter}}');
                AddToLog(out_directory,idx_string,['\\subsection*{\\underline{' title '}}']);
            end
            
            
            fig_index_updated = fig_index + 1;
            idx_string        = ['idx_' num2str(fig_index_updated,'%03i')];
            AddToLog(out_directory,idx_string,subtitle,fig_name);
            %AddToLog('Run_Output/',idx_string,'LarssonVsSourbron','Larsson_vs_Sourbron.png');


end