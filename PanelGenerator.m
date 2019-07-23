classdef PanelGenerator
methods(Static)
    function confusion
        res = load('confusion_mat_one_session.mat');
        nt = res.num_trials/2;
        C_pct = res.C/nt*100;
        C_s_pct = res.C_s/nt*100;
        figure;
        imagesc(C_pct, [0 100]);
        xlabel 'Predicted bin'
        ylabel 'Correct bin'
        axis equal;
        nb = size(res.C,1);
        xlim([1 nb] + [-0.5 0.5]);
        set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        line([nb nb]/2+0.5, ylim, 'Color', 'w');
        line(xlim, [nb nb]/2+0.5, 'Color', 'w');
        h_ = colorbar; ylabel(h_, 'Confusion (%)', 'Rotation', 270);
        Utils.specific_format('confusion');
        Utils.printto('figure1_pdf', 'confusion_unshuffled');
        
        figure;
        imagesc(C_s_pct, [0 100]);
        xlabel 'Predicted bin'
        ylabel 'Correct bin'
        axis equal;
        nb = size(res.C,1);
        xlim([1 nb] + [-0.5 0.5]);
        set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        line([nb nb]/2+0.5, ylim, 'Color', 'w');
        line(xlim, [nb nb]/2+0.5, 'Color', 'w');
        h_ = colorbar; ylabel(h_, 'Confusion (%)', 'Rotation', 270);
        Utils.specific_format('confusion');
        Utils.printto('figure1_pdf', 'confusion_shuffled');
        
        figure;
        imagesc(C_s_pct - C_pct);
        xlabel 'Predicted bin'
        ylabel 'Correct bin'
        axis equal;
        nb = size(res.C,1);
        xlim([1 nb] + [-0.5 0.5]);
        set(gca, 'XTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        set(gca, 'YTickLabel', []);%{'60 cm', '120 cm', '60 cm', '120 cm'});
        line([nb nb]/2+0.5, ylim, 'Color', 'k');
        line(xlim, [nb nb]/2+0.5, 'Color', 'k');
        colormap(bluewhitered);
        h_ = colorbar; ylabel(h_, 'Confusion difference (%)', 'Rotation', 270);
        Utils.specific_format('confusion');
        Utils.printto('figure1_pdf', 'confusion_diff');
    end
end
end