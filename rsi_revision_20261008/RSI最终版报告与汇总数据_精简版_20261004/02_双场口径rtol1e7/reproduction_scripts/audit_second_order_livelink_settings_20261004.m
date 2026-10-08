function audit_second_order_livelink_settings_20261004
% Read-only saved-model acceptance; no runAll, study.run, mesh.run or save.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'revision_comsol','second_order_nr_rtol1e7_dual_origin_20261004');
addpath('D:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
import com.comsol.model.util.*
mphstart('localhost',str2double(getenv('SECOND_ORDER_COMSOL_PORT')), ...
    getenv('COMSOL_REV_USER'),getenv('COMSOL_REV_PASSWORD'));
records=struct([]);
for origin={'entrance','midpoint'}
    for name={'H0_control','H0p5','H1','H2'}
        d=fullfile(out,origin{1},name{1}); model=mphload(fullfile(d,'settings_model.mph'));
        cleanup=onCleanup(@() ModelUtil.remove(model.tag()));
        t=model.sol('sol1').feature('t1'); s=model.study('std1').feature('time');
        row=struct('field_reference',string(origin{1}),'case_name',string(name{1}));
        for property={'rtol','timemethod','tstepsgenalpha','maxstepconstraintgenalpha', ...
                'maxstepgenalpha','rhoinf','atolglobal','atolglobalmethod', ...
                'atolglobalvaluemethod','atolglobalfactor','ewtrescale'}
            row.(property{1})=string(t.getString(property{1}));
        end
        row.study_usertol=string(s.getString('usertol'));
        row.study_rtol=string(s.getString('rtol'));
        row.relativistic_correction=string(model.component('comp1').physics('cpt'). ...
            prop('RelativisticCorrection').getString('RelativisticCorrection'));
        row.new_solves=0;
        assert(str2double(row.rtol)==1e-7 && str2double(row.study_rtol)==1e-7 && row.study_usertol=="on");
        assert(row.tstepsgenalpha=="free" && str2double(row.maxstepgenalpha)==4e-15 && ...
            row.relativistic_correction=="0" && row.atolglobalvaluemethod=="factor" && ...
            str2double(row.atolglobalfactor)==0.1 && row.atolglobalmethod=="scaled");
        records=[records;row]; %#ok<AGROW>
        writetable(struct2table(records),fullfile(out,'saved_model_livelink_settings_audit.csv'));
        fprintf('READ_ONLY_MODEL_ACCEPTANCE %s/%s PASS\n',origin{1},name{1});
        clear cleanup model;
    end
end
fid=fopen(fullfile(out,'LIVELINK_MODEL_READBACK_COMPLETE.txt'),'w');
fprintf(fid,'Eight saved models read back; zero new solves, zero saves.\n'); fclose(fid);
fprintf('LIVELINK_SAVED_MODELS_ALL_PASS\n');
end
