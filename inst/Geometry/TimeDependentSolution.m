classdef TimeDependentSolution < Solution
   
    properties
      
    end
    
    methods
 
       function obj = TimeDependentSolution(asb, d)
        obj@Solution(asb, d)
        obj.domain = asb.domain;
        obj.id = asb.id_matrix;
        obj.d = d;
        [s1, ~] = size(obj.id);
        cpoints = cell(1,s1);
            for i=1:s1
                cpoints{i} = [d(i,:) 1];
            end
        cpoints = reshape(cpoints, size(obj.points));
        obj.cpoints = cpoints;
       end
       
       function recordSolution(obj, name,framerate)
            figure(1)
            vidFile = VideoWriter(name, 'MPEG-4');
            vidFile.FrameRate = framerate;
            open(vidFile)
            for i = 1:size(obj.d,2)
                gcf = obj.plot_solution(i);
                caxis([0 1])
                xlim([0 1])
                ylim([0 1])
                view(0,90)
                drawnow
                writeVideo(vidFile, F(i));
            end
            close(vidFile)
       end
        
       function snapSolution(obj, name, framerate)
           
           figure(1)
           for i = 1:framerate:size(obj.d, 2)
               gct = obj.plot_solution(i);
               caxis([0 1])
               xlim([0 1])
               ylim([0 1])
               view(0,90)
               drawnow
               saveas(gcf,[name,num2str(i),'.png'])
           end
       end
       
    end
end